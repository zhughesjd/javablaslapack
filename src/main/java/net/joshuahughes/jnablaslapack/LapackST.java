package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackST extends Library
{

	public static LapackST instance = (LapackST) Native.loadLibrary("liblapack",LapackST.class);

/**
*> \brief \b STBCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STBCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            INFO, KD, LDAB, N<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AB( LDAB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STBCON estimates the reciprocal of the condition number of a<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stbcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STBRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STBRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,<br>
*                          LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AB( LDAB, * ), B( LDB, * ), BERR( * ),<br>
*      $                   FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STBRFS provides error bounds and backward error estimates for the<br>
*> solution to a system of linear equations with a triangular band<br>
*> coefficient matrix.<br>
*><br>
*> The solution matrix X must be computed by STBTRS or some other<br>
*> means before entering this routine.  STBRFS does not do iterative<br>
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
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stbrfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STBTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STBTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stbtrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stbtrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stbtrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,<br>
*                          LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STBTRS solves a triangular system of the form<br>
*><br>
*>    A * X = B  or  A**T * X = B,<br>
*><br>
*> where A is a triangular band matrix of order N, and B is an<br>
*> N-by NRHS matrix.  A check is made to verify that A is nonsingular.<br>
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
*>          Specifies the form the system of equations:<br>
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stbtrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b STFSM solves a matrix equation (one operand is a triangular matrix in RFP format).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STFSM + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfsm.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfsm.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfsm.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A,<br>
*                         B, LDB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, DIAG, SIDE, TRANS, UPLO<br>
*       INTEGER            LDB, M, N<br>
*       REAL               ALPHA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( 0: * ), B( 0: LDB-1, 0: * )<br>
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
*> STFSM  solves the matrix equation<br>
*><br>
*>    op( A )*X = alpha*B  or  X*op( A ) = alpha*B<br>
*><br>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T.<br>
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
*>          = 'T':  The Transpose Form of RFP A is stored.<br>
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
*>              TRANS  = 'T' or 't'   op( A ) = A'.<br>
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
*>          ALPHA is REAL<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (NT)<br>
*>           NT = N*(N+1)/2. On entry, the matrix A in RFP Format.<br>
*>           RFP Format is described by TRANSR, UPLO and N as follows:<br>
*>           If TRANSR='N' then RFP A is (0:N,0:K-1) when N is even;<br>
*>           K=N/2. RFP A is (0:N-1,0:K) when N is odd; K=N/2. If<br>
*>           TRANSR = 'T' then RFP is the transpose of RFP A as<br>
*>           defined when TRANSR = 'N'. The contents of RFP A are defined<br>
*>           by UPLO as follows: If UPLO = 'U' the RFP A contains the NT<br>
*>           elements of upper packed A either in normal or<br>
*>           transpose Format. If UPLO = 'L' the RFP A contains<br>
*>           the NT elements of lower packed A either in normal or<br>
*>           transpose Format. The LDA of RFP A is (N+1)/2 when<br>
*>           TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is<br>
*>           even and is N when is odd.<br>
*>           See the Note below for more details. Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, DIMENSION (LDB,N)<br>
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
	public void stfsm_(CHARACTER TRANSR,CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER M,INTEGER N,REAL ALPHA,float[] A,float[] B,INTEGER LDB);
/**
*> \brief \b STFTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STFTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stftri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stftri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stftri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STFTRI( TRANSR, UPLO, DIAG, N, A, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO, DIAG<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STFTRI computes the inverse of a triangular matrix A stored in RFP<br>
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
*>          = 'T':  The Transpose TRANSR of RFP A is stored.<br>
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
*>          A is REAL array, dimension (NT);<br>
*>          NT=N*(N+1)/2. On entry, the triangular factor of a Hermitian<br>
*>          Positive Definite matrix A in RFP format. RFP format is<br>
*>          described by TRANSR, UPLO, and N as follows: If TRANSR = 'N'<br>
*>          then RFP A is (0:N,0:k-1) when N is even; k=N/2. RFP A is<br>
*>          (0:N-1,0:k) when N is odd; k=N/2. IF TRANSR = 'T' then RFP is<br>
*>          the transpose of RFP A as defined when<br>
*>          TRANSR = 'N'. The contents of RFP A are defined by UPLO as<br>
*>          follows: If UPLO = 'U' the RFP A contains the nt elements of<br>
*>          upper packed A; If UPLO = 'L' the RFP A contains the nt<br>
*>          elements of lower packed A. The LDA of RFP A is (N+1)/2 when<br>
*>          TRANSR = 'T'. When TRANSR is 'N' the LDA is N+1 when N is<br>
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
*><br>
*  =====================================================================<br>
*/
	public void stftri_(CHARACTER TRANSR,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] A,INTEGER INFO);
/**
*> \brief \b STFTTP copies a triangular matrix from the rectangular full packed format (TF) to the standard packed format (TP).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STFTTP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfttp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfttp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfttp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STFTTP( TRANSR, UPLO, N, ARF, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( 0: * ), ARF( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STFTTP copies a triangular matrix A from rectangular full packed<br>
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
*>          = 'T':  ARF is in Transpose format;<br>
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
*>          ARF is REAL array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A stored in<br>
*>          RFP format. For a further discussion see Notes below.<br>
*> \endverbatim<br>
*><br>
*> \param[out] AP<br>
*> \verbatim<br>
*>          AP is REAL array, dimension ( N*(N+1)/2 ),<br>
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
*><br>
*  =====================================================================<br>
*/
	public void stfttp_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] ARF,float[] AP,INTEGER INFO);
/**
*> \brief \b STFTTR copies a triangular matrix from the rectangular full packed format (TF) to the standard full format (TR).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STFTTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stfttr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stfttr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stfttr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STFTTR( TRANSR, UPLO, N, ARF, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( 0: LDA-1, 0: * ), ARF( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STFTTR copies a triangular matrix A from rectangular full packed<br>
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
*>          = 'T':  ARF is in Transpose format.<br>
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
*>          The order of the matrices ARF and A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ARF<br>
*> \verbatim<br>
*>          ARF is REAL array, dimension (N*(N+1)/2).<br>
*>          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')<br>
*>          matrix A in RFP format. See the "Notes" below for more<br>
*>          details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
	public void stfttr_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] ARF,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b STGEVC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGEVC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgevc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgevc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgevc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,<br>
*                          LDVL, VR, LDVR, MM, M, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, SIDE<br>
*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       REAL               P( LDP, * ), S( LDS, * ), VL( LDVL, * ),<br>
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
*> STGEVC computes some or all of the right and/or left eigenvectors of<br>
*> a pair of real matrices (S,P), where S is a quasi-triangular matrix<br>
*> and P is upper triangular.  Matrix pairs of this type are produced by<br>
*> the generalized Schur factorization of a matrix pair (A,B):<br>
*><br>
*>    A = Q*S*Z**T,  B = Q*P*Z**T<br>
*><br>
*> as computed by SGGHRD + SHGEQZ.<br>
*><br>
*> The right eigenvector x and the left eigenvector y of (S,P)<br>
*> corresponding to an eigenvalue w are defined by:<br>
*> <br>
*>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,<br>
*> <br>
*> where y**H denotes the conjugate tranpose of y.<br>
*> The eigenvalues are not input to this routine, but are computed<br>
*> directly from the diagonal blocks of S and P.<br>
*> <br>
*> This routine returns the matrices X and/or Y of right and left<br>
*> eigenvectors of (S,P), or the products Z*X and/or Q*Y,<br>
*> where Z and Q are input matrices.<br>
*> If Q and Z are the orthogonal factors from the generalized Schur<br>
*> factorization of a matrix pair (A,B), then Z*X and Q*Y<br>
*> are the matrices of right and left eigenvectors of (A,B).<br>
*> <br>
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
*>          computed.  If w(j) is a real eigenvalue, the corresponding<br>
*>          real eigenvector is computed if SELECT(j) is .TRUE..<br>
*>          If w(j) and w(j+1) are the real and imaginary parts of a<br>
*>          complex eigenvalue, the corresponding complex eigenvector<br>
*>          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,<br>
*>          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is<br>
*>          set to .FALSE..<br>
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
*>          S is REAL array, dimension (LDS,N)<br>
*>          The upper quasi-triangular matrix S from a generalized Schur<br>
*>          factorization, as computed by SHGEQZ.<br>
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
*>          P is REAL array, dimension (LDP,N)<br>
*>          The upper triangular matrix P from a generalized Schur<br>
*>          factorization, as computed by SHGEQZ.<br>
*>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks<br>
*>          of S must be in positive diagonal form.<br>
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
*>          VL is REAL array, dimension (LDVL,MM)<br>
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must<br>
*>          contain an N-by-N matrix Q (usually the orthogonal matrix Q<br>
*>          of left Schur vectors returned by SHGEQZ).<br>
*>          On exit, if SIDE = 'L' or 'B', VL contains:<br>
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);<br>
*>          if HOWMNY = 'B', the matrix Q*Y;<br>
*>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by<br>
*>                      SELECT, stored consecutively in the columns of<br>
*>                      VL, in the same order as their eigenvalues.<br>
*><br>
*>          A complex eigenvector corresponding to a complex eigenvalue<br>
*>          is stored in two consecutive columns, the first holding the<br>
*>          real part, and the second the imaginary part.<br>
*><br>
*>          Not referenced if SIDE = 'R'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of array VL.  LDVL >= 1, and if<br>
*>          SIDE = 'L' or 'B', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is REAL array, dimension (LDVR,MM)<br>
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must<br>
*>          contain an N-by-N matrix Z (usually the orthogonal matrix Z<br>
*>          of right Schur vectors returned by SHGEQZ).<br>
*><br>
*>          On exit, if SIDE = 'R' or 'B', VR contains:<br>
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);<br>
*>          if HOWMNY = 'B' or 'b', the matrix Z*X;<br>
*>          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)<br>
*>                      specified by SELECT, stored consecutively in the<br>
*>                      columns of VR, in the same order as their<br>
*>                      eigenvalues.<br>
*><br>
*>          A complex eigenvector corresponding to a complex eigenvalue<br>
*>          is stored in two consecutive columns, the first holding the<br>
*>          real part and the second the imaginary part.<br>
*>          <br>
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
*>          is set to N.  Each selected real eigenvector occupies one<br>
*>          column and each selected complex eigenvector occupies two<br>
*>          columns.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (6*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex<br>
*>                eigenvalue.<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Allocation of workspace:<br>
*>  ---------- -- ---------<br>
*><br>
*>     WORK( j ) = 1-norm of j-th column of A, above the diagonal<br>
*>     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal<br>
*>     WORK( 2*N+1:3*N ) = real part of eigenvector<br>
*>     WORK( 3*N+1:4*N ) = imaginary part of eigenvector<br>
*>     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector<br>
*>     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector<br>
*><br>
*>  Rowwise vs. columnwise solution methods:<br>
*>  ------- --  ---------- -------- -------<br>
*><br>
*>  Finding a generalized eigenvector consists basically of solving the<br>
*>  singular triangular system<br>
*><br>
*>   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)<br>
*><br>
*>  Consider finding the i-th right eigenvector (assume all eigenvalues<br>
*>  are real). The equation to be solved is:<br>
*>       n                   i<br>
*>  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1<br>
*>      k=j                 k=j<br>
*><br>
*>  where  C = (A - w B)  (The components v(i+1:n) are 0.)<br>
*><br>
*>  The "rowwise" method is:<br>
*><br>
*>  (1)  v(i) := 1<br>
*>  for j = i-1,. . .,1:<br>
*>                          i<br>
*>      (2) compute  s = - sum C(j,k) v(k)   and<br>
*>                        k=j+1<br>
*><br>
*>      (3) v(j) := s / C(j,j)<br>
*><br>
*>  Step 2 is sometimes called the "dot product" step, since it is an<br>
*>  inner product between the j-th row and the portion of the eigenvector<br>
*>  that has been computed so far.<br>
*><br>
*>  The "columnwise" method consists basically in doing the sums<br>
*>  for all the rows in parallel.  As each v(j) is computed, the<br>
*>  contribution of v(j) times the j-th column of C is added to the<br>
*>  partial sums.  Since FORTRAN arrays are stored columnwise, this has<br>
*>  the advantage that at each step, the elements of C that are accessed<br>
*>  are adjacent to one another, whereas with the rowwise method, the<br>
*>  elements accessed at a step are spaced LDS (and LDP) words apart.<br>
*><br>
*>  When finding left eigenvectors, the matrix in question is the<br>
*>  transpose of the one in storage, so the rowwise method then<br>
*>  actually accesses columns of A and B at each step, and so is the<br>
*>  preferred method.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgevc_(CHARACTER SIDE,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,float[] S,INTEGER LDS,float[] P,INTEGER LDP,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,float[] WORK,INTEGER INFO);
/**
*> \brief \b STGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an orthogonal equivalence transformation.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGEX2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgex2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgex2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgex2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,<br>
*                          LDZ, J1, N1, N2, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ, WANTZ<br>
*       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, LWORK, N, N1, N2<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STGEX2 swaps adjacent diagonal blocks (A11, B11) and (A22, B22)<br>
*> of size 1-by-1 or 2-by-2 in an upper (quasi) triangular matrix pair<br>
*> (A, B) by an orthogonal equivalence transformation.<br>
*><br>
*> (A, B) must be in generalized real Schur canonical form (as returned<br>
*> by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2<br>
*> diagonal blocks. B is upper triangular.<br>
*><br>
*> Optionally, the matrices Q and Z of generalized Schur vectors are<br>
*> updated.<br>
*><br>
*>        Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T<br>
*>        Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T<br>
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
*>          A is REAL arrays, dimensions (LDA,N)<br>
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
*>          B is REAL arrays, dimensions (LDB,N)<br>
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
*>          Q is REAL array, dimension (LDZ,N)<br>
*>          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.<br>
*>          On exit, the updated matrix Q.<br>
*>          Not referenced if WANTQ = .FALSE..<br>
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
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          On entry, if WANTZ =.TRUE., the orthogonal matrix Z.<br>
*>          On exit, the updated matrix Z.<br>
*>          Not referenced if WANTZ = .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z. LDZ >= 1.<br>
*>          If WANTZ = .TRUE., LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] J1<br>
*> \verbatim<br>
*>          J1 is INTEGER<br>
*>          The index to the first block (A11, B11). 1 <= J1 <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*>          The order of the first block (A11, B11). N1 = 0, 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N2<br>
*> \verbatim<br>
*>          N2 is INTEGER<br>
*>          The order of the second block (A22, B22). N2 = 0, 1 or 2.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK)).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          LWORK >=  MAX( N*(N2+N1), (N2+N1)*(N2+N1)*2 )<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>            =0: Successful exit<br>
*>            >0: If INFO = 1, the transformed matrix (A, B) would be<br>
*>                too far from generalized Schur form; the blocks are<br>
*>                not swapped and (A, B) and (Q, Z) are unchanged.<br>
*>                The problem of swapping is too ill-conditioned.<br>
*>            <0: If INFO = -16: LWORK is too small. Appropriate value<br>
*>                for LWORK is returned in WORK(1).<br>
*> \endverbatim<br>
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
*> \ingroup realGEauxiliary<br>
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
*> \verbatim<br>
*><br>
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the<br>
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in<br>
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and<br>
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.<br>
*><br>
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified<br>
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition<br>
*>      Estimation: Theory, Algorithms and Software,<br>
*>      Report UMINF - 94.04, Department of Computing Science, Umea<br>
*>      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working<br>
*>      Note 87. To appear in Numerical Algorithms, 1996.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgex2_(LOGICAL WANTQ,LOGICAL WANTZ,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,INTEGER J1,INTEGER N1,INTEGER N2,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b STGEXC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGEXC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgexc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgexc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgexc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,<br>
*                          LDZ, IFST, ILST, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ, WANTZ<br>
*       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STGEXC reorders the generalized real Schur decomposition of a real<br>
*> matrix pair (A,B) using an orthogonal equivalence transformation<br>
*><br>
*>                (A, B) = Q * (A, B) * Z**T,<br>
*><br>
*> so that the diagonal block of (A, B) with row index IFST is moved<br>
*> to row ILST.<br>
*><br>
*> (A, B) must be in generalized real Schur canonical form (as returned<br>
*> by SGGES), i.e. A is block upper triangular with 1-by-1 and 2-by-2<br>
*> diagonal blocks. B is upper triangular.<br>
*><br>
*> Optionally, the matrices Q and Z of generalized Schur vectors are<br>
*> updated.<br>
*><br>
*>        Q(in) * A(in) * Z(in)**T = Q(out) * A(out) * Z(out)**T<br>
*>        Q(in) * B(in) * Z(in)**T = Q(out) * B(out) * Z(out)**T<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the matrix A in generalized real Schur canonical<br>
*>          form.<br>
*>          On exit, the updated matrix A, again in generalized<br>
*>          real Schur canonical form.<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
*>          On entry, the matrix B in generalized real Schur canonical<br>
*>          form (A,B).<br>
*>          On exit, the updated matrix B, again in generalized<br>
*>          real Schur canonical form (A,B).<br>
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
*>          Q is REAL array, dimension (LDZ,N)<br>
*>          On entry, if WANTQ = .TRUE., the orthogonal matrix Q.<br>
*>          On exit, the updated matrix Q.<br>
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
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          On entry, if WANTZ = .TRUE., the orthogonal matrix Z.<br>
*>          On exit, the updated matrix Z.<br>
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
*> \param[in,out] IFST<br>
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
*>          On exit, if IFST pointed on entry to the second row of<br>
*>          a 2-by-2 block, it is changed to point to the first row;<br>
*>          ILST always points to the first row of the block in its<br>
*>          final position (which may differ from its input value by<br>
*>          +1 or -1). 1 <= IFST, ILST <= N.<br>
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
*>          LWORK >= 1 when N <= 1, otherwise LWORK >= 4*N + 16.<br>
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
*>           =0:  successful exit.<br>
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
*> \ingroup realGEcomputational<br>
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
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgexc_(LOGICAL WANTQ,LOGICAL WANTZ,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,INTEGER IFST,INTEGER ILST,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b STGSEN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGSEN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsen.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsen.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsen.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,<br>
*                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL,<br>
*                          PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       LOGICAL            WANTQ, WANTZ<br>
*       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,<br>
*      $                   M, N<br>
*       REAL               PL, PR<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),<br>
*      $                   B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STGSEN reorders the generalized real Schur decomposition of a real<br>
*> matrix pair (A, B) (in terms of an orthonormal equivalence trans-<br>
*> formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues<br>
*> appears in the leading diagonal blocks of the upper quasi-triangular<br>
*> matrix A and the upper triangular B. The leading columns of Q and<br>
*> Z form orthonormal bases of the corresponding left and right eigen-<br>
*> spaces (deflating subspaces). (A, B) must be in generalized real<br>
*> Schur canonical form (as returned by SGGES), i.e. A is block upper<br>
*> triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper<br>
*> triangular.<br>
*><br>
*> STGSEN also computes the generalized eigenvalues<br>
*><br>
*>             w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)<br>
*><br>
*> of the reordered matrix pair (A, B).<br>
*><br>
*> Optionally, STGSEN computes the estimates of reciprocal condition<br>
*> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),<br>
*> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)<br>
*> between the matrix pairs (A11, B11) and (A22,B22) that correspond to<br>
*> the selected cluster and the eigenvalues outside the cluster, resp.,<br>
*> and norms of "projections" onto left and right eigenspaces w.r.t.<br>
*> the selected cluster in the (1,1)-block.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
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
*>          SELECT specifies the eigenvalues in the selected cluster.<br>
*>          To select a real eigenvalue w(j), SELECT(j) must be set to<br>
*>          .TRUE.. To select a complex conjugate pair of eigenvalues<br>
*>          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,<br>
*>          either SELECT(j) or SELECT(j+1) or both must be set to<br>
*>          .TRUE.; a complex conjugate pair of eigenvalues must be<br>
*>          either both included in the cluster or both excluded.<br>
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
*>          A is REAL array, dimension(LDA,N)<br>
*>          On entry, the upper quasi-triangular matrix A, with (A, B) in<br>
*>          generalized real Schur canonical form.<br>
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
*>          B is REAL array, dimension(LDB,N)<br>
*>          On entry, the upper triangular matrix B, with (A, B) in<br>
*>          generalized real Schur canonical form.<br>
*>          On exit, B is overwritten by the reordered matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B. LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is REAL array, dimension (N)<br>
*><br>
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will<br>
*>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i<br>
*>          and BETA(j),j=1,...,N  are the diagonals of the complex Schur<br>
*>          form (S,T) that would result if the 2-by-2 diagonal blocks of<br>
*>          the real generalized Schur form of (A,B) were further reduced<br>
*>          to triangular form using complex unitary transformations.<br>
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if<br>
*>          positive, then the j-th and (j+1)-st eigenvalues are a<br>
*>          complex conjugate pair, with ALPHAI(j+1) negative.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
*>          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.<br>
*>          On exit, Q has been postmultiplied by the left orthogonal<br>
*>          transformation matrix which reorder (A, B); The leading M<br>
*>          columns of Q form orthonormal bases for the specified pair of<br>
*>          left eigenspaces (deflating subspaces).<br>
*>          If WANTQ = .FALSE., Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= 1;<br>
*>          and if WANTQ = .TRUE., LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.<br>
*>          On exit, Z has been postmultiplied by the left orthogonal<br>
*>          transformation matrix which reorder (A, B); The leading M<br>
*>          columns of Z form orthonormal bases for the specified pair of<br>
*>          left eigenspaces (deflating subspaces).<br>
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
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The dimension of the specified pair of left and right eigen-<br>
*>          spaces (deflating subspaces). 0 <= M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] PL<br>
*> \verbatim<br>
*>          PL is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[out] PR<br>
*> \verbatim<br>
*>          PR is REAL<br>
*><br>
*>          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the<br>
*>          reciprocal of the norm of "projections" onto left and right<br>
*>          eigenspaces with respect to the selected cluster.<br>
*>          0 < PL, PR <= 1.<br>
*>          If M = 0 or M = N, PL = PR  = 1.<br>
*>          If IJOB = 0, 2 or 3, PL and PR are not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIF<br>
*> \verbatim<br>
*>          DIF is REAL array, dimension (2).<br>
*>          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.<br>
*>          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on<br>
*>          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based<br>
*>          estimates of Difu and Difl.<br>
*>          If M = 0 or N, DIF(1:2) = F-norm([A, B]).<br>
*>          If IJOB = 0 or 1, DIF is not referenced.<br>
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
*>          The dimension of the array WORK. LWORK >=  4*N+16.<br>
*>          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)).<br>
*>          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)).<br>
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
*>          If IJOB = 1, 2 or 4, LIWORK >=  N+6.<br>
*>          If IJOB = 3 or 5, LIWORK >= MAX(2*M*(N-M), N+6).<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  STGSEN first collects the selected eigenvalues by computing<br>
*>  orthogonal U and W that move them to the top left corner of (A, B).<br>
*>  In other words, the selected eigenvalues are the eigenvalues of<br>
*>  (A11, B11) in:<br>
*><br>
*>              U**T*(A, B)*W = (A11 A12) (B11 B12) n1<br>
*>                              ( 0  A22),( 0  B22) n2<br>
*>                                n1  n2    n1  n2<br>
*><br>
*>  where N = n1+n2 and U**T means the transpose of U. The first n1 columns<br>
*>  of U and W span the specified pair of left and right eigenspaces<br>
*>  (deflating subspaces) of (A, B).<br>
*><br>
*>  If (A, B) has been obtained from the generalized real Schur<br>
*>  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**T, then the<br>
*>  reordered generalized real Schur form of (C, D) is given by<br>
*><br>
*>           (C, D) = (Q*U)*(U**T*(A, B)*W)*(Z*W)**T,<br>
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
*>       Zu = [ kron(In2, A11)  -kron(A22**T, In1) ]<br>
*>            [ kron(In2, B11)  -kron(B22**T, In1) ].<br>
*><br>
*>  Here, Inx is the identity matrix of size nx and A22**T is the<br>
*>  transpose of A22. kron(X, Y) is the Kronecker product between<br>
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
*>  based estimate DIF is not wanted (see SLATDF), then the parameter<br>
*>  IDIFJB (see below) should be changed from 3 to 4 (routine SLATDF<br>
*>  (IJOB = 2 will be used)). See STGSYL for more details.<br>
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
*>      Estimation: Theory, Algorithms and Software,<br>
*>      Report UMINF - 94.04, Department of Computing Science, Umea<br>
*>      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working<br>
*>      Note 87. To appear in Numerical Algorithms, 1996.<br>
*><br>
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working<br>
*>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,<br>
*>      1996.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgsen_(INTEGER IJOB,LOGICAL WANTQ,LOGICAL WANTZ,boolean[] SELECT,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,INTEGER M,REAL PL,REAL PR,float[] DIF,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b STGSJA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGSJA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsja.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsja.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsja.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,<br>
*                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,<br>
*                          Q, LDQ, WORK, NCYCLE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBQ, JOBU, JOBV<br>
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,<br>
*      $                   NCYCLE, P<br>
*       REAL               TOLA, TOLB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), ALPHA( * ), B( LDB, * ),<br>
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
*> STGSJA computes the generalized singular value decomposition (GSVD)<br>
*> of two real upper triangular (or trapezoidal) matrices A and B.<br>
*><br>
*> On entry, it is assumed that matrices A and B have the following<br>
*> forms, which may be obtained by the preprocessing subroutine SGGSVP<br>
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
*>        U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ),<br>
*><br>
*> where U, V and Q are orthogonal matrices.<br>
*> R is a nonsingular upper triangular matrix, and D1 and D2 are<br>
*> ``diagonal'' matrices, which are of the following structures:<br>
*><br>
*> If M-K-L >= 0,<br>
*><br>
*>                     K  L<br>
*>        D1 =     K ( I  0 )<br>
*>                 L ( 0  C )<br>
*>             M-K-L ( 0  0 )<br>
*><br>
*>                   K  L<br>
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
*> The computation of the orthogonal transformation matrices U, V or Q<br>
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
*>          = 'U':  U must contain an orthogonal matrix U1 on entry, and<br>
*>                  the product U1*U is returned;<br>
*>          = 'I':  U is initialized to the unit matrix, and the<br>
*>                  orthogonal matrix U is returned;<br>
*>          = 'N':  U is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBV<br>
*> \verbatim<br>
*>          JOBV is CHARACTER*1<br>
*>          = 'V':  V must contain an orthogonal matrix V1 on entry, and<br>
*>                  the product V1*V is returned;<br>
*>          = 'I':  V is initialized to the unit matrix, and the<br>
*>                  orthogonal matrix V is returned;<br>
*>          = 'N':  V is not computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBQ<br>
*> \verbatim<br>
*>          JOBQ is CHARACTER*1<br>
*>          = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and<br>
*>                  the product Q1*Q is returned;<br>
*>          = 'I':  Q is initialized to the unit matrix, and the<br>
*>                  orthogonal matrix Q is returned;<br>
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
*>          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N)<br>
*>          of A and B, whose GSVD is going to be computed by STGSJA.<br>
*>          See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
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
*>          TOLA is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] TOLB<br>
*> \verbatim<br>
*>          TOLB is REAL<br>
*><br>
*>          TOLA and TOLB are the convergence criteria for the Jacobi-<br>
*>          Kogbetliantz iteration procedure. Generally, they are the<br>
*>          same as used in the preprocessing step, say<br>
*>              TOLA = max(M,N)*norm(A)*MACHEPS,<br>
*>              TOLB = max(P,N)*norm(B)*MACHEPS.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is REAL array, dimension (N)<br>
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
*>          U is REAL array, dimension (LDU,M)<br>
*>          On entry, if JOBU = 'U', U must contain a matrix U1 (usually<br>
*>          the orthogonal matrix returned by SGGSVP).<br>
*>          On exit,<br>
*>          if JOBU = 'I', U contains the orthogonal matrix U;<br>
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
*>          V is REAL array, dimension (LDV,P)<br>
*>          On entry, if JOBV = 'V', V must contain a matrix V1 (usually<br>
*>          the orthogonal matrix returned by SGGSVP).<br>
*>          On exit,<br>
*>          if JOBV = 'I', V contains the orthogonal matrix V;<br>
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
*>          Q is REAL array, dimension (LDQ,N)<br>
*>          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually<br>
*>          the orthogonal matrix returned by SGGSVP).<br>
*>          On exit,<br>
*>          if JOBQ = 'I', Q contains the orthogonal matrix Q;<br>
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
*>          WORK is REAL array, dimension (2*N)<br>
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
*><br>
*> \verbatim<br>
*>  Internal Parameters<br>
*>  ===================<br>
*><br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  STGSJA essentially uses a variant of Kogbetliantz algorithm to reduce<br>
*>  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L<br>
*>  matrix B13 to the form:<br>
*><br>
*>           U1**T *A13*Q1 = C1*R1; V1**T *B13*Q1 = S1*R1,<br>
*><br>
*>  where U1, V1 and Q1 are orthogonal matrix, and Z**T is the transpose<br>
*>  of Z.  C1 and S1 are diagonal matrices satisfying<br>
*><br>
*>                C1**2 + S1**2 = I,<br>
*><br>
*>  and R1 is an L-by-L nonsingular upper triangular matrix.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgsja_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,INTEGER K,INTEGER L,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL TOLA,REAL TOLB,float[] ALPHA,float[] BETA,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] Q,INTEGER LDQ,float[] WORK,INTEGER NCYCLE,INTEGER INFO);
/**
*> \brief \b STGSNA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGSNA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsna.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsna.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsna.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,<br>
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
*       REAL               A( LDA, * ), B( LDB, * ), DIF( * ), S( * ),<br>
*      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STGSNA estimates reciprocal condition numbers for specified<br>
*> eigenvalues and/or eigenvectors of a matrix pair (A, B) in<br>
*> generalized real Schur canonical form (or of any matrix pair<br>
*> (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where<br>
*> Z**T denotes the transpose of Z.<br>
*><br>
*> (A, B) must be in generalized real Schur form (as returned by SGGES),<br>
*> i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal<br>
*> blocks. B is upper triangular.<br>
*><br>
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
*>          for the eigenpair corresponding to a real eigenvalue w(j),<br>
*>          SELECT(j) must be set to .TRUE.. To select condition numbers<br>
*>          corresponding to a complex conjugate pair of eigenvalues w(j)<br>
*>          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be<br>
*>          set to .TRUE..<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The upper quasi-triangular matrix A in the pair (A,B).<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
*>          The upper triangular matrix B in the pair (A,B).<br>
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
*>          VL is REAL array, dimension (LDVL,M)<br>
*>          If JOB = 'E' or 'B', VL must contain left eigenvectors of<br>
*>          (A, B), corresponding to the eigenpairs specified by HOWMNY<br>
*>          and SELECT. The eigenvectors must be stored in consecutive<br>
*>          columns of VL, as returned by STGEVC.<br>
*>          If JOB = 'V', VL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL. LDVL >= 1.<br>
*>          If JOB = 'E' or 'B', LDVL >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VR<br>
*> \verbatim<br>
*>          VR is REAL array, dimension (LDVR,M)<br>
*>          If JOB = 'E' or 'B', VR must contain right eigenvectors of<br>
*>          (A, B), corresponding to the eigenpairs specified by HOWMNY<br>
*>          and SELECT. The eigenvectors must be stored in consecutive<br>
*>          columns ov VR, as returned by STGEVC.<br>
*>          If JOB = 'V', VR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR. LDVR >= 1.<br>
*>          If JOB = 'E' or 'B', LDVR >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (MM)<br>
*>          If JOB = 'E' or 'B', the reciprocal condition numbers of the<br>
*>          selected eigenvalues, stored in consecutive elements of the<br>
*>          array. For a complex conjugate pair of eigenvalues two<br>
*>          consecutive elements of S are set to the same value. Thus<br>
*>          S(j), DIF(j), and the j-th columns of VL and VR all<br>
*>          correspond to the same eigenpair (but not in general the<br>
*>          j-th eigenpair, unless all eigenpairs are selected).<br>
*>          If JOB = 'V', S is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] DIF<br>
*> \verbatim<br>
*>          DIF is REAL array, dimension (MM)<br>
*>          If JOB = 'V' or 'B', the estimated reciprocal condition<br>
*>          numbers of the selected eigenvectors, stored in consecutive<br>
*>          elements of the array. For a complex eigenvector two<br>
*>          consecutive elements of DIF are set to the same value. If<br>
*>          the eigenvalues cannot be reordered to compute DIF(j), DIF(j)<br>
*>          is set to 0; this can only occur when the true value would be<br>
*>          very small anyway.<br>
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
*>          the specified condition numbers; for each selected real<br>
*>          eigenvalue one element is used, and for each selected complex<br>
*>          conjugate pair of eigenvalues, two elements are used.<br>
*>          If HOWMNY = 'A', M is set to N.<br>
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
*>          If JOB = 'V' or 'B' LWORK >= 2*N*(N+2)+16.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N + 6)<br>
*>          If JOB = 'E', IWORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          =0: Successful exit<br>
*>          <0: If INFO = -i, the i-th argument had an illegal value<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The reciprocal of the condition number of a generalized eigenvalue<br>
*>  w = (a, b) is defined as<br>
*><br>
*>       S(w) = (|u**TAv|**2 + |u**TBv|**2)**(1/2) / (norm(u)*norm(v))<br>
*><br>
*>  where u and v are the left and right eigenvectors of (A, B)<br>
*>  corresponding to w; |z| denotes the absolute value of the complex<br>
*>  number, and norm(u) denotes the 2-norm of the vector u.<br>
*>  The pair (a, b) corresponds to an eigenvalue w = a/b (= u**TAv/u**TBv)<br>
*>  of the matrix pair (A, B). If both a and b equal zero, then (A B) is<br>
*>  singular and S(I) = -1 is returned.<br>
*><br>
*>  An approximate error bound on the chordal distance between the i-th<br>
*>  computed generalized eigenvalue w and the corresponding exact<br>
*>  eigenvalue lambda is<br>
*><br>
*>       chord(w, lambda) <= EPS * norm(A, B) / S(I)<br>
*><br>
*>  where EPS is the machine precision.<br>
*><br>
*>  The reciprocal of the condition number DIF(i) of right eigenvector u<br>
*>  and left eigenvector v corresponding to the generalized eigenvalue w<br>
*>  is defined as follows:<br>
*><br>
*>  a) If the i-th eigenvalue w = (a,b) is real<br>
*><br>
*>     Suppose U and V are orthogonal transformations such that<br>
*><br>
*>              U**T*(A, B)*V  = (S, T) = ( a   *  ) ( b  *  )  1<br>
*>                                        ( 0  S22 ),( 0 T22 )  n-1<br>
*>                                          1  n-1     1 n-1<br>
*><br>
*>     Then the reciprocal condition number DIF(i) is<br>
*><br>
*>                Difl((a, b), (S22, T22)) = sigma-min( Zl ),<br>
*><br>
*>     where sigma-min(Zl) denotes the smallest singular value of the<br>
*>     2(n-1)-by-2(n-1) matrix<br>
*><br>
*>         Zl = [ kron(a, In-1)  -kron(1, S22) ]<br>
*>              [ kron(b, In-1)  -kron(1, T22) ] .<br>
*><br>
*>     Here In-1 is the identity matrix of size n-1. kron(X, Y) is the<br>
*>     Kronecker product between the matrices X and Y.<br>
*><br>
*>     Note that if the default method for computing DIF(i) is wanted<br>
*>     (see SLATDF), then the parameter DIFDRI (see below) should be<br>
*>     changed from 3 to 4 (routine SLATDF(IJOB = 2 will be used)).<br>
*>     See STGSYL for more details.<br>
*><br>
*>  b) If the i-th and (i+1)-th eigenvalues are complex conjugate pair,<br>
*><br>
*>     Suppose U and V are orthogonal transformations such that<br>
*><br>
*>              U**T*(A, B)*V = (S, T) = ( S11  *   ) ( T11  *  )  2<br>
*>                                       ( 0    S22 ),( 0    T22) n-2<br>
*>                                         2    n-2     2    n-2<br>
*><br>
*>     and (S11, T11) corresponds to the complex conjugate eigenvalue<br>
*>     pair (w, conjg(w)). There exist unitary matrices U1 and V1 such<br>
*>     that<br>
*><br>
*>       U1**T*S11*V1 = ( s11 s12 ) and U1**T*T11*V1 = ( t11 t12 )<br>
*>                      (  0  s22 )                    (  0  t22 )<br>
*><br>
*>     where the generalized eigenvalues w = s11/t11 and<br>
*>     conjg(w) = s22/t22.<br>
*><br>
*>     Then the reciprocal condition number DIF(i) is bounded by<br>
*><br>
*>         min( d1, max( 1, |real(s11)/real(s22)| )*d2 )<br>
*><br>
*>     where, d1 = Difl((s11, t11), (s22, t22)) = sigma-min(Z1), where<br>
*>     Z1 is the complex 2-by-2 matrix<br>
*><br>
*>              Z1 =  [ s11  -s22 ]<br>
*>                    [ t11  -t22 ],<br>
*><br>
*>     This is done by computing (using real arithmetic) the<br>
*>     roots of the characteristical polynomial det(Z1**T * Z1 - lambda I),<br>
*>     where Z1**T denotes the transpose of Z1 and det(X) denotes<br>
*>     the determinant of X.<br>
*><br>
*>     and d2 is an upper bound on Difl((S11, T11), (S22, T22)), i.e. an<br>
*>     upper bound on sigma-min(Z2), where Z2 is (2n-2)-by-(2n-2)<br>
*><br>
*>              Z2 = [ kron(S11**T, In-2)  -kron(I2, S22) ]<br>
*>                   [ kron(T11**T, In-2)  -kron(I2, T22) ]<br>
*><br>
*>     Note that if the default method for computing DIF is wanted (see<br>
*>     SLATDF), then the parameter DIFDRI (see below) should be changed<br>
*>     from 3 to 4 (routine SLATDF(IJOB = 2 will be used)). See STGSYL<br>
*>     for more details.<br>
*><br>
*>  For each eigenvalue/vector specified by SELECT, DIF stores a<br>
*>  Frobenius norm-based estimate of Difl.<br>
*><br>
*>  An approximate error bound for the i-th computed eigenvector VL(i) or<br>
*>  VR(i) is given by<br>
*><br>
*>             EPS * norm(A, B) / DIF(i).<br>
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
*>      Estimation: Theory, Algorithms and Software,<br>
*>      Report UMINF - 94.04, Department of Computing Science, Umea<br>
*>      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working<br>
*>      Note 87. To appear in Numerical Algorithms, 1996.<br>
*><br>
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working<br>
*>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,<br>
*>      No 1, 1996.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgsna_(CHARACTER JOB,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] S,float[] DIF,INTEGER MM,INTEGER M,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STGSY2 solves the generalized Sylvester equation (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGSY2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsy2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsy2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsy2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,<br>
*                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,<br>
*                          IWORK, PQ, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N,<br>
*      $                   PQ<br>
*       REAL               RDSCAL, RDSUM, SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ),<br>
*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STGSY2 solves the generalized Sylvester equation:<br>
*><br>
*>             A * R - L * B = scale * C                (1)<br>
*>             D * R - L * E = scale * F,<br>
*><br>
*> using Level 1 and 2 BLAS. where R and L are unknown M-by-N matrices,<br>
*> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,<br>
*> N-by-N and M-by-N, respectively, with real entries. (A, D) and (B, E)<br>
*> must be in generalized Schur canonical form, i.e. A, B are upper<br>
*> quasi triangular and D, E are upper triangular. The solution (R, L)<br>
*> overwrites (C, F). 0 <= SCALE <= 1 is an output scaling factor<br>
*> chosen to avoid overflow.<br>
*><br>
*> In matrix notation solving equation (1) corresponds to solve<br>
*> Z*x = scale*b, where Z is defined as<br>
*><br>
*>        Z = [ kron(In, A)  -kron(B**T, Im) ]             (2)<br>
*>            [ kron(In, D)  -kron(E**T, Im) ],<br>
*><br>
*> Ik is the identity matrix of size k and X**T is the transpose of X.<br>
*> kron(X, Y) is the Kronecker product between the matrices X and Y.<br>
*> In the process of solving (1), we solve a number of such systems<br>
*> where Dim(In), Dim(In) = 1 or 2.<br>
*><br>
*> If TRANS = 'T', solve the transposed system Z**T*y = scale*b for y,<br>
*> which is equivalent to solve for R and L in<br>
*><br>
*>             A**T * R  + D**T * L   = scale * C           (3)<br>
*>             R  * B**T + L  * E**T  = scale * -F<br>
*><br>
*> This case is used to compute an estimate of Dif[(A, D), (B, E)] =<br>
*> sigma_min(Z) using reverse communicaton with SLACON.<br>
*><br>
*> STGSY2 also (IJOB >= 1) contributes to the computation in STGSYL<br>
*> of an upper bound on the separation between to matrix pairs. Then<br>
*> the input (A, D), (B, E) are sub-pencils of the matrix pair in<br>
*> STGSYL. See STGSYL for details.<br>
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
*>          = 0: solve (1) only.<br>
*>          = 1: A contribution from this subsystem to a Frobenius<br>
*>               norm-based estimate of the separation between two matrix<br>
*>               pairs is computed. (look ahead strategy is used).<br>
*>          = 2: A contribution from this subsystem to a Frobenius<br>
*>               norm-based estimate of the separation between two matrix<br>
*>               pairs is computed. (SGECON on sub-systems is used.)<br>
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
*>          A is REAL array, dimension (LDA, M)<br>
*>          On entry, A contains an upper quasi triangular matrix.<br>
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
*>          B is REAL array, dimension (LDB, N)<br>
*>          On entry, B contains an upper quasi triangular matrix.<br>
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
*>          C is REAL array, dimension (LDC, N)<br>
*>          On entry, C contains the right-hand-side of the first matrix<br>
*>          equation in (1).<br>
*>          On exit, if IJOB = 0, C has been overwritten by the<br>
*>          solution R.<br>
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
*>          D is REAL array, dimension (LDD, M)<br>
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
*>          E is REAL array, dimension (LDE, N)<br>
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
*>          F is REAL array, dimension (LDF, N)<br>
*>          On entry, F contains the right-hand-side of the second matrix<br>
*>          equation in (1).<br>
*>          On exit, if IJOB = 0, F has been overwritten by the<br>
*>          solution L.<br>
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
*>          SCALE is REAL<br>
*>          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions<br>
*>          R and L (C and F on entry) will hold the solutions to a<br>
*>          slightly perturbed system but the input matrices A, B, D and<br>
*>          E have not been changed. If SCALE = 0, R and L will hold the<br>
*>          solutions to the homogeneous system with C = F = 0. Normally,<br>
*>          SCALE = 1.<br>
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
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (M+N+2)<br>
*> \endverbatim<br>
*><br>
*> \param[out] PQ<br>
*> \verbatim<br>
*>          PQ is INTEGER<br>
*>          On exit, the number of subsystems (of size 2-by-2, 4-by-4 and<br>
*>          8-by-8) solved by this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          On exit, if INFO is set to<br>
*>            =0: Successful exit<br>
*>            <0: If INFO = -i, the i-th argument had an illegal value.<br>
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
*> \ingroup realSYauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,<br>
*>     Umea University, S-901 87 Umea, Sweden.<br>
*<br>
*  =====================================================================<br>
*/
	public void stgsy2_(CHARACTER TRANS,INTEGER IJOB,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,INTEGER LDC,float[] D,INTEGER LDD,float[] E,INTEGER LDE,float[] F,INTEGER LDF,REAL SCALE,REAL RDSUM,REAL RDSCAL,int[] IWORK,INTEGER PQ,INTEGER INFO);
/**
*> \brief \b STGSYL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STGSYL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsyl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsyl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsyl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,<br>
*                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF,<br>
*      $                   LWORK, M, N<br>
*       REAL               DIF, SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * ),<br>
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
*> STGSYL solves the generalized Sylvester equation:<br>
*><br>
*>             A * R - L * B = scale * C                 (1)<br>
*>             D * R - L * E = scale * F<br>
*><br>
*> where R and L are unknown m-by-n matrices, (A, D), (B, E) and<br>
*> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,<br>
*> respectively, with real entries. (A, D) and (B, E) must be in<br>
*> generalized (real) Schur canonical form, i.e. A, B are upper quasi<br>
*> triangular and D, E are upper triangular.<br>
*><br>
*> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output<br>
*> scaling factor chosen to avoid overflow.<br>
*><br>
*> In matrix notation (1) is equivalent to solve  Zx = scale b, where<br>
*> Z is defined as<br>
*><br>
*>            Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)<br>
*>                [ kron(In, D)  -kron(E**T, Im) ].<br>
*><br>
*> Here Ik is the identity matrix of size k and X**T is the transpose of<br>
*> X. kron(X, Y) is the Kronecker product between the matrices X and Y.<br>
*><br>
*> If TRANS = 'T', STGSYL solves the transposed system Z**T*y = scale*b,<br>
*> which is equivalent to solve for R and L in<br>
*><br>
*>             A**T * R + D**T * L = scale * C           (3)<br>
*>             R * B**T + L * E**T = scale * -F<br>
*><br>
*> This case (TRANS = 'T') is used to compute an one-norm-based estimate<br>
*> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)<br>
*> and (B,E), using SLACON.<br>
*><br>
*> If IJOB >= 1, STGSYL computes a Frobenius norm-based estimate<br>
*> of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the<br>
*> reciprocal of the smallest singular value of Z. See [1-2] for more<br>
*> information.<br>
*><br>
*> This is a level 3 BLAS algorithm.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N', solve the generalized Sylvester equation (1).<br>
*>          = 'T', solve the 'transposed' system (3).<br>
*> \endverbatim<br>
*><br>
*> \param[in] IJOB<br>
*> \verbatim<br>
*>          IJOB is INTEGER<br>
*>          Specifies what kind of functionality to be performed.<br>
*>           =0: solve (1) only.<br>
*>           =1: The functionality of 0 and 3.<br>
*>           =2: The functionality of 0 and 4.<br>
*>           =3: Only an estimate of Dif[(A,D), (B,E)] is computed.<br>
*>               (look ahead strategy IJOB  = 1 is used).<br>
*>           =4: Only an estimate of Dif[(A,D), (B,E)] is computed.<br>
*>               ( SGECON on sub-systems is used ).<br>
*>          Not referenced if TRANS = 'T'.<br>
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
*>          A is REAL array, dimension (LDA, M)<br>
*>          The upper quasi triangular matrix A.<br>
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
*>          B is REAL array, dimension (LDB, N)<br>
*>          The upper quasi triangular matrix B.<br>
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
*>          C is REAL array, dimension (LDC, N)<br>
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
*>          D is REAL array, dimension (LDD, M)<br>
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
*>          E is REAL array, dimension (LDE, N)<br>
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
*>          F is REAL array, dimension (LDF, N)<br>
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
*>          DIF is REAL<br>
*>          On exit DIF is the reciprocal of a lower bound of the<br>
*>          reciprocal of the Dif-function, i.e. DIF is an upper bound of<br>
*>          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2).<br>
*>          IF IJOB = 0 or TRANS = 'T', DIF is not touched.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCALE<br>
*> \verbatim<br>
*>          SCALE is REAL<br>
*>          On exit SCALE is the scaling factor in (1) or (3).<br>
*>          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,<br>
*>          to a slightly perturbed system but the input matrices A, B, D<br>
*>          and E have not been changed. If SCALE = 0, C and F hold the<br>
*>          solutions R and L, respectively, to the homogeneous system<br>
*>          with C = F = 0. Normally, SCALE = 1.<br>
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
*>          IWORK is INTEGER array, dimension (M+N+6)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>            =0: successful exit<br>
*>            <0: If INFO = -i, the i-th argument had an illegal value.<br>
*>            >0: (A, D) and (B, E) have common or close eigenvalues.<br>
*> \endverbatim<br>
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
*> \ingroup realSYcomputational<br>
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
*>  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software<br>
*>      for Solving the Generalized Sylvester Equation and Estimating the<br>
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,<br>
*>      Department of Computing Science, Umea University, S-901 87 Umea,<br>
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working<br>
*>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,<br>
*>      No 1, 1996.<br>
*><br>
*>  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester<br>
*>      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.<br>
*>      Appl., 15(4):1045-1060, 1994<br>
*><br>
*>  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with<br>
*>      Condition Estimators for Solving the Generalized Sylvester<br>
*>      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,<br>
*>      July 1989, pp 745-751.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stgsyl_(CHARACTER TRANS,INTEGER IJOB,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,INTEGER LDC,float[] D,INTEGER LDD,float[] E,INTEGER LDE,float[] F,INTEGER LDF,REAL SCALE,REAL DIF,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STPCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            INFO, N<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPCON estimates the reciprocal of the condition number of a packed<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stpcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] AP,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STPMQRT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPMQRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpmqrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpmqrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpmqrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,<br>
*                           A, LDA, B, LDB, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER SIDE, TRANS<br>
*       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), <br>
*      $          WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPMQRT applies a real orthogonal matrix Q obtained from a <br>
*> "triangular-pentagonal" real block reflector H to a general<br>
*> real matrix C, which consists of two blocks A and B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'L': apply Q or Q^T from the Left;<br>
*>          = 'R': apply Q or Q^T from the Right.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N':  No transpose, apply Q;<br>
*>          = 'T':  Transpose, apply Q^T.<br>
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
*>          V is REAL array, dimension (LDA,K)<br>
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
*>          T is REAL array, dimension (LDT,K)<br>
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
*>          A is REAL array, dimension<br>
*>          (LDA,N) if SIDE = 'L' or <br>
*>          (LDA,K) if SIDE = 'R'<br>
*>          On entry, the K-by-N or M-by-K matrix A.<br>
*>          On exit, A is overwritten by the corresponding block of <br>
*>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
*>          On entry, the M-by-N matrix B.<br>
*>          On exit, B is overwritten by the corresponding block of<br>
*>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.<br>
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
*>          WORK is REAL array. The dimension of WORK is<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup realOTHERcomputational<br>
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
*>  The real orthogonal matrix Q is formed from V and T.<br>
*><br>
*>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.<br>
*><br>
*>  If TRANS='T' and SIDE='L', C is on exit replaced with Q^T * C.<br>
*><br>
*>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.<br>
*><br>
*>  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q^T.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stpmqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,INTEGER NB,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] WORK,INTEGER INFO);
/**
*> \brief \b STPQRT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPQRT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpqrt.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpqrt.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpqrt.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INFO, LDA, LDB, LDT, N, M, L, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPQRT computes a blocked QR factorization of a real <br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
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
*>          T is REAL array, dimension (LDT,N)<br>
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
*>          WORK is REAL array, dimension (NB*N)<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void stpqrt_(INTEGER M,INTEGER N,INTEGER L,INTEGER NB,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] T,INTEGER LDT,float[] WORK,INTEGER INFO);
/**
*> \brief \b STPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPQRT2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpqrt2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpqrt2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpqrt2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER   INFO, LDA, LDB, LDT, N, M, L<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL   A( LDA, * ), B( LDB, * ), T( LDT, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPQRT2 computes a QR factorization of a real "triangular-pentagonal"<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
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
*>          T is REAL array, dimension (LDT,N)<br>
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
*> \ingroup realOTHERcomputational<br>
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
*>               H = I - W * T * W^H<br>
*><br>
*>  where W^H is the conjugate transpose of W and T is the upper triangular<br>
*>  factor of the block reflector.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stpqrt2_(INTEGER M,INTEGER N,INTEGER L,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] T,INTEGER LDT,INTEGER INFO);
/**
*> \brief \b STPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex matrix, which is composed of two blocks.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPRFB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stprfb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stprfb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stprfb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, <br>
*                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER DIRECT, SIDE, STOREV, TRANS<br>
*       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL   A( LDA, * ), B( LDB, * ), T( LDT, * ), <br>
*      $          V( LDV, * ), WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPRFB applies a real "triangular-pentagonal" block reflector H or its <br>
*> conjugate transpose H^H to a real matrix C, which is composed of two <br>
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
*>          = 'L': apply H or H^H from the Left<br>
*>          = 'R': apply H or H^H from the Right<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>          = 'N': apply H (No transpose)<br>
*>          = 'C': apply H^H (Conjugate transpose)<br>
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
*>          V is REAL array, dimension<br>
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
*>          T is REAL array, dimension (LDT,K)<br>
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
*>          A is REAL array, dimension<br>
*>          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R'<br>
*>          On entry, the K-by-N or M-by-K matrix A.<br>
*>          On exit, A is overwritten by the corresponding block of <br>
*>          H*C or H^H*C or C*H or C*H^H.  See Futher Details.<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
*>          On entry, the M-by-N matrix B.<br>
*>          On exit, B is overwritten by the corresponding block of<br>
*>          H*C or H^H*C or C*H or C*H^H.  See Further Details.<br>
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
*>          WORK is REAL array, dimension<br>
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
*> \ingroup realOTHERauxiliary<br>
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
	public void stprfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] WORK,INTEGER LDWORK);
/**
*> \brief \b STPRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stprfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stprfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stprfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,<br>
*                          FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AP( * ), B( LDB, * ), BERR( * ), FERR( * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPRFS provides error bounds and backward error estimates for the<br>
*> solution to a system of linear equations with a triangular packed<br>
*> coefficient matrix.<br>
*><br>
*> The solution matrix X must be computed by STPTRS or some other<br>
*> means before entering this routine.  STPRFS does not do iterative<br>
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
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          If DIAG = 'U', the diagonal elements of A are not referenced<br>
*>          and are assumed to be 1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stprfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STPTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stptri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stptri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stptri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPTRI( UPLO, DIAG, N, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPTRI computes the inverse of a real upper or lower triangular<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void stptri_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] AP,INTEGER INFO);
/**
*> \brief \b STPTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stptrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stptrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stptrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPTRS solves a triangular system of the form<br>
*><br>
*>    A * X = B  or  A**T * X = B,<br>
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
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangular matrix A, packed columnwise in<br>
*>          a linear array.  The j-th column of A is stored in the array<br>
*>          AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stptrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b STPTTF copies a triangular matrix from the standard packed format (TP) to the rectangular full packed format (TF).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPTTF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpttf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpttf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpttf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPTTF( TRANSR, UPLO, N, AP, ARF, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( 0: * ), ARF( 0: * )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPTTF copies a triangular matrix A from standard packed format (TP)<br>
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
*>          = 'T':  ARF in Conjugate-transpose format is wanted.<br>
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
*>          AP is REAL array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A, packed<br>
*>          columnwise in a linear array. The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ARF<br>
*> \verbatim<br>
*>          ARF is REAL array, dimension ( N*(N+1)/2 ),<br>
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
*><br>
*  =====================================================================<br>
*/
	public void stpttf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] AP,float[] ARF,INTEGER INFO);
/**
*> \brief \b STPTTR copies a triangular matrix from the standard packed format (TP) to the standard full format (TR).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STPTTR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpttr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpttr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpttr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPTTR( UPLO, N, AP, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPTTR copies a triangular matrix A from standard packed format (TP)<br>
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
*>          AP is REAL array, dimension ( N*(N+1)/2 ),<br>
*>          On entry, the upper or lower triangular matrix A, packed<br>
*>          columnwise in a linear array. The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension ( LDA, N )<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void stpttr_(CHARACTER UPLO,INTEGER N,float[] AP,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b STRCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, NORM, UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRCON estimates the reciprocal of the condition number of a<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          RCOND is REAL<br>
*>          The reciprocal of the condition number of the matrix A,<br>
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STREVC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STREVC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strevc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strevc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strevc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,<br>
*                          LDVR, MM, M, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, SIDE<br>
*       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       REAL               T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STREVC computes some or all of the right and/or left eigenvectors of<br>
*> a real upper quasi-triangular matrix T.<br>
*> Matrices of this type are produced by the Schur factorization of<br>
*> a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.<br>
*> <br>
*> The right eigenvector x and the left eigenvector y of T corresponding<br>
*> to an eigenvalue w are defined by:<br>
*> <br>
*>    T*x = w*x,     (y**T)*T = w*(y**T)<br>
*> <br>
*> where y**T denotes the transpose of y.<br>
*> The eigenvalues are not input to this routine, but are read directly<br>
*> from the diagonal blocks of T.<br>
*> <br>
*> This routine returns the matrices X and/or Y of right and left<br>
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an<br>
*> input matrix.  If Q is the orthogonal factor that reduces a matrix<br>
*> A to Schur form T, then Q*X and Q*Y are the matrices of right and<br>
*> left eigenvectors of A.<br>
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
*>                  backtransformed by the matrices in VR and/or VL;<br>
*>          = 'S':  compute selected right and/or left eigenvectors,<br>
*>                  as indicated by the logical array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be<br>
*>          computed.<br>
*>          If w(j) is a real eigenvalue, the corresponding real<br>
*>          eigenvector is computed if SELECT(j) is .TRUE..<br>
*>          If w(j) and w(j+1) are the real and imaginary parts of a<br>
*>          complex eigenvalue, the corresponding complex eigenvector is<br>
*>          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and<br>
*>          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to<br>
*>          .FALSE..<br>
*>          Not referenced if HOWMNY = 'A' or 'B'.<br>
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
*>          T is REAL array, dimension (LDT,N)<br>
*>          The upper quasi-triangular matrix T in Schur canonical form.<br>
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
*>          VL is REAL array, dimension (LDVL,MM)<br>
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must<br>
*>          contain an N-by-N matrix Q (usually the orthogonal matrix Q<br>
*>          of Schur vectors returned by SHSEQR).<br>
*>          On exit, if SIDE = 'L' or 'B', VL contains:<br>
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*Y;<br>
*>          if HOWMNY = 'S', the left eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VL, in the same order as their<br>
*>                           eigenvalues.<br>
*>          A complex eigenvector corresponding to a complex eigenvalue<br>
*>          is stored in two consecutive columns, the first holding the<br>
*>          real part, and the second the imaginary part.<br>
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
*>          VR is REAL array, dimension (LDVR,MM)<br>
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must<br>
*>          contain an N-by-N matrix Q (usually the orthogonal matrix Q<br>
*>          of Schur vectors returned by SHSEQR).<br>
*>          On exit, if SIDE = 'R' or 'B', VR contains:<br>
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*X;<br>
*>          if HOWMNY = 'S', the right eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VR, in the same order as their<br>
*>                           eigenvalues.<br>
*>          A complex eigenvector corresponding to a complex eigenvalue<br>
*>          is stored in two consecutive columns, the first holding the<br>
*>          real part and the second the imaginary part.<br>
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
*>          used to store the eigenvectors.<br>
*>          If HOWMNY = 'A' or 'B', M is set to N.<br>
*>          Each selected real eigenvector occupies one column and each<br>
*>          selected complex eigenvector occupies two columns.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (3*N)<br>
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
	public void strevc_(CHARACTER SIDE,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,float[] T,INTEGER LDT,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,float[] WORK,INTEGER INFO);
/**
*> \brief \b STREVC3<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STREVC3 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strevc3.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strevc3.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strevc3.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL,<br>
*                           VR, LDVR, MM, M, WORK, LWORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, SIDE<br>
*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       REAL   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STREVC3 computes some or all of the right and/or left eigenvectors of<br>
*> a real upper quasi-triangular matrix T.<br>
*> Matrices of this type are produced by the Schur factorization of<br>
*> a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.<br>
*><br>
*> The right eigenvector x and the left eigenvector y of T corresponding<br>
*> to an eigenvalue w are defined by:<br>
*><br>
*>    T*x = w*x,     (y**T)*T = w*(y**T)<br>
*><br>
*> where y**T denotes the transpose of the vector y.<br>
*> The eigenvalues are not input to this routine, but are read directly<br>
*> from the diagonal blocks of T.<br>
*><br>
*> This routine returns the matrices X and/or Y of right and left<br>
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an<br>
*> input matrix. If Q is the orthogonal factor that reduces a matrix<br>
*> A to Schur form T, then Q*X and Q*Y are the matrices of right and<br>
*> left eigenvectors of A.<br>
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
*>                  backtransformed by the matrices in VR and/or VL;<br>
*>          = 'S':  compute selected right and/or left eigenvectors,<br>
*>                  as indicated by the logical array SELECT.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be<br>
*>          computed.<br>
*>          If w(j) is a real eigenvalue, the corresponding real<br>
*>          eigenvector is computed if SELECT(j) is .TRUE..<br>
*>          If w(j) and w(j+1) are the real and imaginary parts of a<br>
*>          complex eigenvalue, the corresponding complex eigenvector is<br>
*>          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and<br>
*>          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to<br>
*>          .FALSE..<br>
*>          Not referenced if HOWMNY = 'A' or 'B'.<br>
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
*>          T is REAL array, dimension (LDT,N)<br>
*>          The upper quasi-triangular matrix T in Schur canonical form.<br>
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
*>          VL is REAL array, dimension (LDVL,MM)<br>
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must<br>
*>          contain an N-by-N matrix Q (usually the orthogonal matrix Q<br>
*>          of Schur vectors returned by SHSEQR).<br>
*>          On exit, if SIDE = 'L' or 'B', VL contains:<br>
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*Y;<br>
*>          if HOWMNY = 'S', the left eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VL, in the same order as their<br>
*>                           eigenvalues.<br>
*>          A complex eigenvector corresponding to a complex eigenvalue<br>
*>          is stored in two consecutive columns, the first holding the<br>
*>          real part, and the second the imaginary part.<br>
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
*>          VR is REAL array, dimension (LDVR,MM)<br>
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must<br>
*>          contain an N-by-N matrix Q (usually the orthogonal matrix Q<br>
*>          of Schur vectors returned by SHSEQR).<br>
*>          On exit, if SIDE = 'R' or 'B', VR contains:<br>
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;<br>
*>          if HOWMNY = 'B', the matrix Q*X;<br>
*>          if HOWMNY = 'S', the right eigenvectors of T specified by<br>
*>                           SELECT, stored consecutively in the columns<br>
*>                           of VR, in the same order as their<br>
*>                           eigenvalues.<br>
*>          A complex eigenvector corresponding to a complex eigenvalue<br>
*>          is stored in two consecutive columns, the first holding the<br>
*>          real part and the second the imaginary part.<br>
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
*>          Each selected real eigenvector occupies one column and each<br>
*>          selected complex eigenvector occupies two columns.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of array WORK. LWORK >= max(1,3*N).<br>
*>          For optimum performance, LWORK >= N + 2*N*NB, where NB is<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2011<br>
*<br>
*  @generated from dtrevc3.f, fortran d -> s, Tue Apr 19 01:47:44 2016<br>
*<br>
*> \ingroup realOTHERcomputational<br>
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
	public void strevc3_(CHARACTER SIDE,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,float[] T,INTEGER LDT,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b STREXC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STREXC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strexc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strexc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strexc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ<br>
*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N<br>
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
*> STREXC reorders the real Schur factorization of a real matrix<br>
*> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is<br>
*> moved to row ILST.<br>
*><br>
*> The real Schur form T is reordered by an orthogonal similarity<br>
*> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors<br>
*> is updated by postmultiplying it with Z.<br>
*><br>
*> T must be in Schur canonical form (as returned by SHSEQR), that is,<br>
*> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each<br>
*> 2-by-2 diagonal block has its diagonal elements equal and its<br>
*> off-diagonal elements of opposite sign.<br>
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
*>          T is REAL array, dimension (LDT,N)<br>
*>          On entry, the upper quasi-triangular matrix T, in Schur<br>
*>          Schur canonical form.<br>
*>          On exit, the reordered upper quasi-triangular matrix, again<br>
*>          in Schur canonical form.<br>
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
*>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.<br>
*>          On exit, if COMPQ = 'V', Q has been postmultiplied by the<br>
*>          orthogonal transformation matrix Z which reorders T.<br>
*>          If COMPQ = 'N', Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IFST<br>
*> \verbatim<br>
*>          IFST is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] ILST<br>
*> \verbatim<br>
*>          ILST is INTEGER<br>
*><br>
*>          Specify the reordering of the diagonal blocks of T.<br>
*>          The block with row index IFST is moved to row ILST, by a<br>
*>          sequence of transpositions between adjacent blocks.<br>
*>          On exit, if IFST pointed on entry to the second row of a<br>
*>          2-by-2 block, it is changed to point to the first row; ILST<br>
*>          always points to the first row of the block in its final<br>
*>          position (which may differ from its input value by +1 or -1).<br>
*>          1 <= IFST <= N; 1 <= ILST <= N.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          = 1:  two adjacent blocks were too close to swap (the problem<br>
*>                is very ill-conditioned); T may have been partially<br>
*>                reordered, and ILST points to the first row of the<br>
*>                current position of the block being moved.<br>
*> \endverbatim<br>
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
	public void strexc_(CHARACTER COMPQ,INTEGER N,float[] T,INTEGER LDT,float[] Q,INTEGER LDQ,INTEGER IFST,INTEGER ILST,float[] WORK,INTEGER INFO);
/**
*> \brief \b STRRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,<br>
*                          LDX, FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRRFS provides error bounds and backward error estimates for the<br>
*> solution to a system of linear equations with a triangular<br>
*> coefficient matrix.<br>
*><br>
*> The solution matrix X must be computed by STRTRS or some other<br>
*> means before entering this routine.  STRRFS does not do iterative<br>
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
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strrfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STRSEN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRSEN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsen.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsen.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsen.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI,<br>
*                          M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, JOB<br>
*       INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N<br>
*       REAL               S, SEP<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IWORK( * )<br>
*       REAL               Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ),<br>
*      $                   WR( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRSEN reorders the real Schur factorization of a real matrix<br>
*> A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in<br>
*> the leading diagonal blocks of the upper quasi-triangular matrix T,<br>
*> and the leading columns of Q form an orthonormal basis of the<br>
*> corresponding right invariant subspace.<br>
*><br>
*> Optionally the routine computes the reciprocal condition numbers of<br>
*> the cluster of eigenvalues and/or the invariant subspace.<br>
*><br>
*> T must be in Schur canonical form (as returned by SHSEQR), that is,<br>
*> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each<br>
*> 2-by-2 diagonal block has its diagonal elements equal and its<br>
*> off-diagonal elements of opposite sign.<br>
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
*>          select a real eigenvalue w(j), SELECT(j) must be set to<br>
*>          .TRUE.. To select a complex conjugate pair of eigenvalues<br>
*>          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,<br>
*>          either SELECT(j) or SELECT(j+1) or both must be set to<br>
*>          .TRUE.; a complex conjugate pair of eigenvalues must be<br>
*>          either both included in the cluster or both excluded.<br>
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
*>          On exit, T is overwritten by the reordered matrix T, again in<br>
*>          Schur canonical form, with the selected eigenvalues in the<br>
*>          leading diagonal blocks.<br>
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
*>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.<br>
*>          On exit, if COMPQ = 'V', Q has been postmultiplied by the<br>
*>          orthogonal transformation matrix which reorders T; the<br>
*>          leading M columns of Q form an orthonormal basis for the<br>
*>          specified invariant subspace.<br>
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
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is REAL array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is REAL array, dimension (N)<br>
*><br>
*>          The real and imaginary parts, respectively, of the reordered<br>
*>          eigenvalues of T. The eigenvalues are stored in the same<br>
*>          order as on the diagonal of T, with WR(i) = T(i,i) and, if<br>
*>          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and<br>
*>          WI(i+1) = -WI(i). Note that if a complex eigenvalue is<br>
*>          sufficiently ill-conditioned, then its value may differ<br>
*>          significantly from its value before reordering.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The dimension of the specified invariant subspace.<br>
*>          0 < = M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] S<br>
*> \verbatim<br>
*>          S is REAL<br>
*>          If JOB = 'E' or 'B', S is a lower bound on the reciprocal<br>
*>          condition number for the selected cluster of eigenvalues.<br>
*>          S cannot underestimate the true reciprocal condition number<br>
*>          by more than a factor of sqrt(N). If M = 0 or N, S = 1.<br>
*>          If JOB = 'N' or 'V', S is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SEP<br>
*> \verbatim<br>
*>          SEP is REAL<br>
*>          If JOB = 'V' or 'B', SEP is the estimated reciprocal<br>
*>          condition number of the specified invariant subspace. If<br>
*>          M = 0 or N, SEP = norm(T).<br>
*>          If JOB = 'N' or 'E', SEP is not referenced.<br>
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
*>          If JOB = 'N', LWORK >= max(1,N);<br>
*>          if JOB = 'E', LWORK >= max(1,M*(N-M));<br>
*>          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).<br>
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
*>          The dimension of the array IWORK.<br>
*>          If JOB = 'N' or 'E', LIWORK >= 1;<br>
*>          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).<br>
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
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          = 1: reordering of T failed because some eigenvalues are too<br>
*>               close to separate (the problem is very ill-conditioned);<br>
*>               T may have been partially reordered, and WR and WI<br>
*>               contain the eigenvalues in the same order as in T; S and<br>
*>               SEP (if requested) are set to zero.<br>
*> \endverbatim<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  STRSEN first collects the selected eigenvalues by computing an<br>
*>  orthogonal transformation Z to move them to the top left corner of T.<br>
*>  In other words, the selected eigenvalues are the eigenvalues of T11<br>
*>  in:<br>
*><br>
*>          Z**T * T * Z = ( T11 T12 ) n1<br>
*>                         (  0  T22 ) n2<br>
*>                            n1  n2<br>
*><br>
*>  where N = n1+n2 and Z**T means the transpose of Z. The first n1 columns<br>
*>  of Z span the specified invariant subspace of T.<br>
*><br>
*>  If T has been obtained from the real Schur factorization of a matrix<br>
*>  A = Q*T*Q**T, then the reordered real Schur factorization of A is given<br>
*>  by A = (Q*Z)*(Z**T*T*Z)*(Q*Z)**T, and the first n1 columns of Q*Z span<br>
*>  the corresponding invariant subspace of A.<br>
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
	public void strsen_(CHARACTER JOB,CHARACTER COMPQ,boolean[] SELECT,INTEGER N,float[] T,INTEGER LDT,float[] Q,INTEGER LDQ,float[] WR,float[] WI,INTEGER M,REAL S,REAL SEP,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b STRSNA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRSNA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsna.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsna.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsna.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,<br>
*                          LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          HOWMNY, JOB<br>
*       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IWORK( * )<br>
*       REAL               S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ),<br>
*      $                   VR( LDVR, * ), WORK( LDWORK, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRSNA estimates reciprocal condition numbers for specified<br>
*> eigenvalues and/or right eigenvectors of a real upper<br>
*> quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q<br>
*> orthogonal).<br>
*><br>
*> T must be in Schur canonical form (as returned by SHSEQR), that is,<br>
*> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each<br>
*> 2-by-2 diagonal block has its diagonal elements equal and its<br>
*> off-diagonal elements of opposite sign.<br>
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
*>          for the eigenpair corresponding to a real eigenvalue w(j),<br>
*>          SELECT(j) must be set to .TRUE.. To select condition numbers<br>
*>          corresponding to a complex conjugate pair of eigenvalues w(j)<br>
*>          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be<br>
*>          set to .TRUE..<br>
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
*>          T is REAL array, dimension (LDT,N)<br>
*>          The upper quasi-triangular matrix T, in Schur canonical form.<br>
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
*>          VL is REAL array, dimension (LDVL,M)<br>
*>          If JOB = 'E' or 'B', VL must contain left eigenvectors of T<br>
*>          (or of any Q*T*Q**T with Q orthogonal), corresponding to the<br>
*>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors<br>
*>          must be stored in consecutive columns of VL, as returned by<br>
*>          SHSEIN or STREVC.<br>
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
*>          VR is REAL array, dimension (LDVR,M)<br>
*>          If JOB = 'E' or 'B', VR must contain right eigenvectors of T<br>
*>          (or of any Q*T*Q**T with Q orthogonal), corresponding to the<br>
*>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors<br>
*>          must be stored in consecutive columns of VR, as returned by<br>
*>          SHSEIN or STREVC.<br>
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
*>          S is REAL array, dimension (MM)<br>
*>          If JOB = 'E' or 'B', the reciprocal condition numbers of the<br>
*>          selected eigenvalues, stored in consecutive elements of the<br>
*>          array. For a complex conjugate pair of eigenvalues two<br>
*>          consecutive elements of S are set to the same value. Thus<br>
*>          S(j), SEP(j), and the j-th columns of VL and VR all<br>
*>          correspond to the same eigenpair (but not in general the<br>
*>          j-th eigenpair, unless all eigenpairs are selected).<br>
*>          If JOB = 'V', S is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SEP<br>
*> \verbatim<br>
*>          SEP is REAL array, dimension (MM)<br>
*>          If JOB = 'V' or 'B', the estimated reciprocal condition<br>
*>          numbers of the selected eigenvectors, stored in consecutive<br>
*>          elements of the array. For a complex eigenvector two<br>
*>          consecutive elements of SEP are set to the same value. If<br>
*>          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)<br>
*>          is set to 0; this can only occur when the true value would be<br>
*>          very small anyway.<br>
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
*>          WORK is REAL array, dimension (LDWORK,N+6)<br>
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
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (2*(N-1))<br>
*>          If JOB = 'E', IWORK is not referenced.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The reciprocal of the condition number of an eigenvalue lambda is<br>
*>  defined as<br>
*><br>
*>          S(lambda) = |v**T*u| / (norm(u)*norm(v))<br>
*><br>
*>  where u and v are the right and left eigenvectors of T corresponding<br>
*>  to lambda; v**T denotes the transpose of v, and norm(u)<br>
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
	public void strsna_(CHARACTER JOB,CHARACTER HOWMNY,boolean[] SELECT,INTEGER N,float[] T,INTEGER LDT,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] S,float[] SEP,INTEGER MM,INTEGER M,float[] WORK,INTEGER LDWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b STRSYL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRSYL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strsyl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strsyl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strsyl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,<br>
*                          LDC, SCALE, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANA, TRANB<br>
*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N<br>
*       REAL               SCALE<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRSYL solves the real Sylvester matrix equation:<br>
*><br>
*>    op(A)*X + X*op(B) = scale*C or<br>
*>    op(A)*X - X*op(B) = scale*C,<br>
*><br>
*> where op(A) = A or A**T, and  A and B are both upper quasi-<br>
*> triangular. A is M-by-M and B is N-by-N; the right hand side C and<br>
*> the solution X are M-by-N; and scale is an output scale factor, set<br>
*> <= 1 to avoid overflow in X.<br>
*><br>
*> A and B must be in Schur canonical form (as returned by SHSEQR), that<br>
*> is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;<br>
*> each 2-by-2 diagonal block has its diagonal elements equal and its<br>
*> off-diagonal elements of opposite sign.<br>
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
*>          = 'T': op(A) = A**T (Transpose)<br>
*>          = 'C': op(A) = A**H (Conjugate transpose = Transpose)<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANB<br>
*> \verbatim<br>
*>          TRANB is CHARACTER*1<br>
*>          Specifies the option op(B):<br>
*>          = 'N': op(B) = B    (No transpose)<br>
*>          = 'T': op(B) = B**T (Transpose)<br>
*>          = 'C': op(B) = B**H (Conjugate transpose = Transpose)<br>
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
*>          A is REAL array, dimension (LDA,M)<br>
*>          The upper quasi-triangular matrix A, in Schur canonical form.<br>
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
*>          B is REAL array, dimension (LDB,N)<br>
*>          The upper quasi-triangular matrix B, in Schur canonical form.<br>
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
*>          C is REAL array, dimension (LDC,N)<br>
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
*>          SCALE is REAL<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strsyl_(CHARACTER TRANA,CHARACTER TRANB,INTEGER ISGN,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,INTEGER LDC,REAL SCALE,INTEGER INFO);
/**
*> \brief \b STRTI2 computes the inverse of a triangular matrix (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRTI2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strti2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strti2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strti2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRTI2( UPLO, DIAG, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, UPLO<br>
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
*> STRTI2 computes the inverse of a real upper or lower triangular<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strti2_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b STRTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strtri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strtri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strtri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRTRI( UPLO, DIAG, N, A, LDA, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, UPLO<br>
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
*> STRTRI computes the inverse of a real upper or lower triangular<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strtri_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
/**
*> \brief \b STRTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strtrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strtrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strtrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG, TRANS, UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
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
*> STRTRS solves a triangular system of the form<br>
*><br>
*>    A * X = B  or  A**T * X = B,<br>
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
*>          = 'N':  A * X = B  (No transpose)<br>
*>          = 'T':  A**T * X = B  (Transpose)<br>
*>          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strtrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b STRTTF copies a triangular matrix from the standard full format (TR) to the rectangular full packed format (TF).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRTTF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strttf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strttf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strttf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRTTF( TRANSR, UPLO, N, A, LDA, ARF, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANSR, UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( 0: LDA-1, 0: * ), ARF( 0: * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRTTF copies a triangular matrix A from standard full format (TR)<br>
*> to rectangular full packed format (TF) .<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSR<br>
*> \verbatim<br>
*>          TRANSR is CHARACTER*1<br>
*>          = 'N':  ARF in Normal form is wanted;<br>
*>          = 'T':  ARF in Transpose form is wanted.<br>
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
*>          The order of the matrix A. N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N).<br>
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
*>          The leading dimension of the matrix A. LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ARF<br>
*> \verbatim<br>
*>          ARF is REAL array, dimension (NT).<br>
*>          NT=N*(N+1)/2. On exit, the triangular matrix A in RFP format.<br>
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
	public void strttf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] ARF,INTEGER INFO);
/**
*> \brief \b STRTTP copies a triangular matrix from the standard full format (TR) to the standard packed format (TP).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STRTTP + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strttp.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strttp.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strttp.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRTTP( UPLO, N, A, LDA, AP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRTTP copies a triangular matrix A from full format (TR) to standard<br>
*> packed format (TP).<br>
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
*>          The order of the matrices AP and A.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \param[out] AP<br>
*> \verbatim<br>
*>          AP is REAL array, dimension (N*(N+1)/2<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void strttp_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AP,INTEGER INFO);
/**
*> \brief \b STZRZF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download STZRZF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stzrzf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stzrzf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stzrzf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, LWORK, M, N<br>
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
*> STZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A<br>
*> to upper triangular form by means of orthogonal transformations.<br>
*><br>
*> The upper trapezoidal matrix A is factored as<br>
*><br>
*>    A = ( R  0 ) * Z,<br>
*><br>
*> where Z is an N-by-N orthogonal matrix and R is an M-by-M upper<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the leading M-by-N upper trapezoidal part of the<br>
*>          array A must contain the matrix to be factorized.<br>
*>          On exit, the leading M-by-M upper triangular part of A<br>
*>          contains the upper triangular matrix R, and elements M+1 to<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*>  The N-by-N matrix Z can be computed by<br>
*><br>
*>     Z =  Z(1)*Z(2)* ... *Z(M)<br>
*><br>
*>  where each N-by-N Z(k) is given by<br>
*><br>
*>     Z(k) = I - tau(k)*v(k)*v(k)**T<br>
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
	public void stzrzf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);

}