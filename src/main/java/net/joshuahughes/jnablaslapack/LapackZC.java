package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZC extends Library
{

	public static LapackZC instance = (LapackZC) Native.loadLibrary("liblapack",LapackZC.class);

/**
*> \brief <b> ZCGESV computes the solution to system of linear equations A * X = B for GE matrices</b> (mixed precision with iterative refinement)<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZCGESV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zcgesv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zcgesv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zcgesv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZCGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK,<br>
*                          SWORK, RWORK, ITER, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX            SWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( N, * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZCGESV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.<br>
*><br>
*> ZCGESV first attempts to factorize the matrix in COMPLEX and use this<br>
*> factorization within an iterative refinement procedure to produce a<br>
*> solution with COMPLEX*16 normwise backward error quality (see below).<br>
*> If the approach fails the method switches to a COMPLEX*16<br>
*> factorization and solve.<br>
*><br>
*> The iterative refinement is not going to be a winning strategy if<br>
*> the ratio COMPLEX performance over COMPLEX*16 performance is too<br>
*> small. A reasonable strategy should take the number of right-hand<br>
*> sides and the size of the matrix into account. This might be done<br>
*> with a call to ILAENV in the future. Up to now, we always try<br>
*> iterative refinement.<br>
*><br>
*> The iterative refinement process is stopped if<br>
*>     ITER > ITERMAX<br>
*> or for all the RHS we have:<br>
*>     RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX<br>
*> where<br>
*>     o ITER is the number of the current iteration in the iterative<br>
*>       refinement process<br>
*>     o RNRM is the infinity-norm of the residual<br>
*>     o XNRM is the infinity-norm of the solution<br>
*>     o ANRM is the infinity-operator-norm of the matrix A<br>
*>     o EPS is the machine epsilon returned by DLAMCH('Epsilon')<br>
*> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00<br>
*> respectively.<br>
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
*>          A is COMPLEX*16 array,<br>
*>          dimension (LDA,N)<br>
*>          On entry, the N-by-N coefficient matrix A.<br>
*>          On exit, if iterative refinement has been successfully used<br>
*>          (INFO.EQ.0 and ITER.GE.0, see description below), then A is<br>
*>          unchanged, if double precision factorization has been used<br>
*>          (INFO.EQ.0 and ITER.LT.0, see description below), then the<br>
*>          array A contains the factors L and U from the factorization<br>
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
*>          Corresponds either to the single precision factorization<br>
*>          (if INFO.EQ.0 and ITER.GE.0) or the double precision<br>
*>          factorization (if INFO.EQ.0 and ITER.LT.0).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          If INFO = 0, the N-by-NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N*NRHS)<br>
*>          This array is used to hold the residual vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SWORK<br>
*> \verbatim<br>
*>          SWORK is COMPLEX array, dimension (N*(N+NRHS))<br>
*>          This array is used to use the single precision matrix and the<br>
*>          right-hand sides or solutions in single precision.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ITER<br>
*> \verbatim<br>
*>          ITER is INTEGER<br>
*>          < 0: iterative refinement has failed, COMPLEX*16<br>
*>               factorization has been performed<br>
*>               -1 : the routine fell back to full precision for<br>
*>                    implementation- or machine-specific reasons<br>
*>               -2 : narrowing the precision induced an overflow,<br>
*>                    the routine fell back to full precision<br>
*>               -3 : failure of CGETRF<br>
*>               -31: stop the iterative refinement after the 30th<br>
*>                    iterations<br>
*>          > 0: iterative refinement has been successfully used.<br>
*>               Returns the number of iterations<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, U(i,i) computed in COMPLEX*16 is exactly<br>
*>                zero.  The factorization has been completed, but the<br>
*>                factor U is exactly singular, so the solution<br>
*>                could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16GEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zcgesv_(INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] WORK,float[] SWORK,double[] RWORK,INTEGER ITER,INTEGER INFO);
/**
*> \brief <b> ZCPOSV computes the solution to system of linear equations A * X = B for PO matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZCPOSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zcposv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zcposv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zcposv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZCPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK,<br>
*                          SWORK, RWORK, ITER, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX            SWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( N, * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZCPOSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian positive definite matrix and X and B<br>
*> are N-by-NRHS matrices.<br>
*><br>
*> ZCPOSV first attempts to factorize the matrix in COMPLEX and use this<br>
*> factorization within an iterative refinement procedure to produce a<br>
*> solution with COMPLEX*16 normwise backward error quality (see below).<br>
*> If the approach fails the method switches to a COMPLEX*16<br>
*> factorization and solve.<br>
*><br>
*> The iterative refinement is not going to be a winning strategy if<br>
*> the ratio COMPLEX performance over COMPLEX*16 performance is too<br>
*> small. A reasonable strategy should take the number of right-hand<br>
*> sides and the size of the matrix into account. This might be done<br>
*> with a call to ILAENV in the future. Up to now, we always try<br>
*> iterative refinement.<br>
*><br>
*> The iterative refinement process is stopped if<br>
*>     ITER > ITERMAX<br>
*> or for all the RHS we have:<br>
*>     RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX<br>
*> where<br>
*>     o ITER is the number of the current iteration in the iterative<br>
*>       refinement process<br>
*>     o RNRM is the infinity-norm of the residual<br>
*>     o XNRM is the infinity-norm of the solution<br>
*>     o ANRM is the infinity-operator-norm of the matrix A<br>
*>     o EPS is the machine epsilon returned by DLAMCH('Epsilon')<br>
*> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00<br>
*> respectively.<br>
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
*>          A is COMPLEX*16 array,<br>
*>          dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A. If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          Note that the imaginary parts of the diagonal<br>
*>          elements need not be set and are assumed to be zero.<br>
*><br>
*>          On exit, if iterative refinement has been successfully used<br>
*>          (INFO.EQ.0 and ITER.GE.0, see description below), then A is<br>
*>          unchanged, if double precision factorization has been used<br>
*>          (INFO.EQ.0 and ITER.LT.0, see description below), then the<br>
*>          array A contains the factor U or L from the Cholesky<br>
*>          factorization A = U**H*U or A = L*L**H.<br>
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
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          If INFO = 0, the N-by-NRHS solution matrix X.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.  LDX >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N*NRHS)<br>
*>          This array is used to hold the residual vectors.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SWORK<br>
*> \verbatim<br>
*>          SWORK is COMPLEX array, dimension (N*(N+NRHS))<br>
*>          This array is used to use the single precision matrix and the<br>
*>          right-hand sides or solutions in single precision.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ITER<br>
*> \verbatim<br>
*>          ITER is INTEGER<br>
*>          < 0: iterative refinement has failed, COMPLEX*16<br>
*>               factorization has been performed<br>
*>               -1 : the routine fell back to full precision for<br>
*>                    implementation- or machine-specific reasons<br>
*>               -2 : narrowing the precision induced an overflow,<br>
*>                    the routine fell back to full precision<br>
*>               -3 : failure of CPOTRF<br>
*>               -31: stop the iterative refinement after the 30th<br>
*>                    iterations<br>
*>          > 0: iterative refinement has been successfully used.<br>
*>               Returns the number of iterations<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the leading minor of order i of<br>
*>                (COMPLEX*16) A is not positive definite, so the<br>
*>                factorization could not be completed, and the solution<br>
*>                has not been computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16POsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zcposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] WORK,float[] SWORK,double[] RWORK,INTEGER ITER,INTEGER INFO);

}