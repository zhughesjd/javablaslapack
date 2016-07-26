package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZS extends Library
{

	public static LapackZS instance = (LapackZS) Native.loadLibrary("liblapack",LapackZS.class);

/**
*> \brief \b ZSPCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPCON estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a complex symmetric packed matrix A using the<br>
*> factorization A = U*D*U**T or A = L*D*L**T computed by ZSPTRF.<br>
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
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSPTRF, stored as a<br>
*>          packed triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSPTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is DOUBLE PRECISION<br>
*>          The 1-norm of the original matrix A.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
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
	public void zspcon_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSPMV computes a matrix-vector product for complex vectors using a complex symmetric packed matrix<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspmv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspmv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspmv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INCX, INCY, N<br>
*       COMPLEX*16         ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPMV  performs the matrix-vector operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n symmetric matrix, supplied in packed form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the matrix A is supplied in the packed<br>
*>           array AP as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   The upper triangular part of A is<br>
*>                                  supplied in AP.<br>
*><br>
*>              UPLO = 'L' or 'l'   The lower triangular part of A is<br>
*>                                  supplied in AP.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension at least<br>
*>           ( ( N*( N + 1 ) )/2 ).<br>
*>           Before entry, with UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on.<br>
*>           Before entry, with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the N-<br>
*>           element vector x.<br>
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
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y. On exit, Y is overwritten by the updated<br>
*>           vector y.<br>
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
*> \ingroup complex16OTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void zspmv_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] AP,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZSPR performs the symmetrical rank-1 update of a complex symmetric packed matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPR( UPLO, N, ALPHA, X, INCX, AP )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INCX, N<br>
*       COMPLEX*16         ALPHA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPR    performs the symmetric rank 1 operation<br>
*><br>
*>    A := alpha*x*x**H + A,<br>
*><br>
*> where alpha is a complex scalar, x is an n element vector and A is an<br>
*> n by n symmetric matrix, supplied in packed form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the matrix A is supplied in the packed<br>
*>           array AP as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   The upper triangular part of A is<br>
*>                                  supplied in AP.<br>
*><br>
*>              UPLO = 'L' or 'l'   The lower triangular part of A is<br>
*>                                  supplied in AP.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the N-<br>
*>           element vector x.<br>
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
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension at least<br>
*>           ( ( N*( N + 1 ) )/2 ).<br>
*>           Before entry, with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry, with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the lower triangular part of the<br>
*>           updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set, they are assumed to be zero, and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
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
*  =====================================================================<br>
*/
	public void zspr_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] AP);
/**
*> \brief \b ZSPRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsprfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsprfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsprfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,<br>
*                          FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         AFP( * ), AP( * ), B( LDB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is symmetric indefinite<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangle of the symmetric matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFP<br>
*> \verbatim<br>
*>          AFP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The factored form of the matrix A.  AFP contains the block<br>
*>          diagonal matrix D and the multipliers used to obtain the<br>
*>          factor U or L from the factorization A = U*D*U**T or<br>
*>          A = L*D*L**T as computed by ZSPTRF, stored as a packed<br>
*>          triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSPTRF.<br>
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
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by ZSPTRS.<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZSPSV computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N symmetric matrix stored in packed format and X<br>
*> and B are N-by-NRHS matrices.<br>
*><br>
*> The diagonal pivoting method is used to factor A as<br>
*>    A = U * D * U**T,  if UPLO = 'U', or<br>
*>    A = L * D * L**T,  if UPLO = 'L',<br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, D is symmetric and block diagonal with 1-by-1<br>
*> and 2-by-2 diagonal blocks.  The factored form of A is then used to<br>
*> solve the system of equations A * X = B.<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          See below for further details.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by ZSPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D, as<br>
*>          determined by ZSPTRF.  If IPIV(k) > 0, then rows and columns<br>
*>          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1<br>
*>          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,<br>
*>          then rows and columns k-1 and -IPIV(k) were interchanged and<br>
*>          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and<br>
*>          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and<br>
*>          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2<br>
*>          diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization<br>
*>                has been completed, but the block diagonal matrix D is<br>
*>                exactly singular, so the solution could not be<br>
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
*> \ingroup complex16OTHERsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The packed storage scheme is illustrated by the following example<br>
*>  when N = 4, UPLO = 'U':<br>
*><br>
*>  Two-dimensional storage of the symmetric matrix A:<br>
*><br>
*>     a11 a12 a13 a14<br>
*>         a22 a23 a24<br>
*>             a33 a34     (aij = aji)<br>
*>                 a44<br>
*><br>
*>  Packed storage of the upper triangle of A:<br>
*><br>
*>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zspsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> ZSPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X,<br>
*                          LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          FACT, UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         AFP( * ), AP( * ), B( LDB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPSVX uses the diagonal pivoting factorization A = U*D*U**T or<br>
*> A = L*D*L**T to compute the solution to a complex system of linear<br>
*> equations A * X = B, where A is an N-by-N symmetric matrix stored<br>
*> in packed format and X and B are N-by-NRHS matrices.<br>
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
*> 1. If FACT = 'N', the diagonal pivoting method is used to factor A as<br>
*>       A = U * D * U**T,  if UPLO = 'U', or<br>
*>       A = L * D * L**T,  if UPLO = 'L',<br>
*>    where U (or L) is a product of permutation and unit upper (lower)<br>
*>    triangular matrices and D is symmetric and block diagonal with<br>
*>    1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*> 2. If some D(i,i)=0, so that D is exactly singular, then the routine<br>
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
*>          = 'F':  On entry, AFP and IPIV contain the factored form<br>
*>                  of A.  AP, AFP and IPIV will not be modified.<br>
*>          = 'N':  The matrix A will be copied to AFP and factored.<br>
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
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangle of the symmetric matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          See below for further details.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AFP<br>
*> \verbatim<br>
*>          AFP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          If FACT = 'F', then AFP is an input argument and on entry<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by ZSPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*><br>
*>          If FACT = 'N', then AFP is an output argument and on exit<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by ZSPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          If FACT = 'F', then IPIV is an input argument and on entry<br>
*>          contains details of the interchanges and the block structure<br>
*>          of D, as determined by ZSPTRF.<br>
*>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>          interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and<br>
*>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =<br>
*>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were<br>
*>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*><br>
*>          If FACT = 'N', then IPIV is an output argument and on exit<br>
*>          contains details of the interchanges and the block structure<br>
*>          of D, as determined by ZSPTRF.<br>
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
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is<br>
*>                <= N:  D(i,i) is exactly zero.  The factorization<br>
*>                       has been completed but the factor D is exactly<br>
*>                       singular, so the solution and error bounds could<br>
*>                       not be computed. RCOND = 0 is returned.<br>
*>                = N+1: D is nonsingular, but RCOND is less than machine<br>
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
*> \ingroup complex16OTHERsolve<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The packed storage scheme is illustrated by the following example<br>
*>  when N = 4, UPLO = 'U':<br>
*><br>
*>  Two-dimensional storage of the symmetric matrix A:<br>
*><br>
*>     a11 a12 a13 a14<br>
*>         a22 a23 a24<br>
*>             a33 a34     (aij = aji)<br>
*>                 a44<br>
*><br>
*>  Packed storage of the upper triangle of A:<br>
*><br>
*>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zspsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZSPTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsptrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsptrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsptrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPTRF( UPLO, N, AP, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPTRF computes the factorization of a complex symmetric matrix A<br>
*> stored in packed format using the Bunch-Kaufman diagonal pivoting<br>
*> method:<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is symmetric and block diagonal with<br>
*> 1-by-1 and 2-by-2 diagonal blocks.<br>
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
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L, stored as a packed triangular<br>
*>          matrix overwriting A (see below for further details).<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D.<br>
*>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>          interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and<br>
*>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =<br>
*>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were<br>
*>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular, and division by zero will occur if it<br>
*>               is used to solve a system of equations.<br>
*> \endverbatim<br>
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
*>  5-96 - Based on modifications by J. Lewis, Boeing Computer Services<br>
*>         Company<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**T, where<br>
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,<br>
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to<br>
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    v    0   )   k-s<br>
*>     U(k) =  (   0    I    0   )   s<br>
*>             (   0    0    I   )   n-k<br>
*>                k-s   s   n-k<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).<br>
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),<br>
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).<br>
*><br>
*>  If UPLO = 'L', then A = L*D*L**T, where<br>
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,<br>
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to<br>
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    0     0   )  k-1<br>
*>     L(k) =  (   0    I     0   )  s<br>
*>             (   0    v     I   )  n-k-s+1<br>
*>                k-1   s  n-k-s+1<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).<br>
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),<br>
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zsptrf_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,INTEGER INFO);
/**
*> \brief \b ZSPTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsptri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsptri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsptri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPTRI( UPLO, N, AP, IPIV, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPTRI computes the inverse of a complex symmetric indefinite matrix<br>
*> A in packed storage using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by ZSPTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          On entry, the block diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by ZSPTRF,<br>
*>          stored as a packed triangular matrix.<br>
*><br>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original<br>
*>          matrix, stored as a packed triangular matrix. The j-th column<br>
*>          of inv(A) is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L',<br>
*>             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSPTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its<br>
*>               inverse could not be computed.<br>
*> \endverbatim<br>
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
	public void zsptri_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSPTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSPTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsptrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsptrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsptrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSPTRS solves a system of linear equations A*X = B with a complex<br>
*> symmetric matrix A stored in packed format using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by ZSPTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSPTRF, stored as a<br>
*>          packed triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSPTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZSTEDC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSTEDC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstedc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstedc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstedc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK,<br>
*                          LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPZ<br>
*       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   D( * ), E( * ), RWORK( * )<br>
*       COMPLEX*16         WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSTEDC computes all eigenvalues and, optionally, eigenvectors of a<br>
*> symmetric tridiagonal matrix using the divide and conquer method.<br>
*> The eigenvectors of a full or band complex Hermitian matrix can also<br>
*> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this<br>
*> matrix to tridiagonal form.<br>
*><br>
*> This code makes very mild assumptions about floating point<br>
*> arithmetic. It will work on machines with a guard digit in<br>
*> add/subtract, or on those binary machines without guard digits<br>
*> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.<br>
*> It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.  See DLAED3 for details.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only.<br>
*>          = 'I':  Compute eigenvectors of tridiagonal matrix also.<br>
*>          = 'V':  Compute eigenvectors of original Hermitian matrix<br>
*>                  also.  On entry, Z contains the unitary matrix used<br>
*>                  to reduce the original matrix to tridiagonal form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The dimension of the symmetric tridiagonal matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the diagonal elements of the tridiagonal matrix.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, the subdiagonal elements of the tridiagonal matrix.<br>
*>          On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ,N)<br>
*>          On entry, if COMPZ = 'V', then Z contains the unitary<br>
*>          matrix used in the reduction to tridiagonal form.<br>
*>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the<br>
*>          orthonormal eigenvectors of the original Hermitian matrix,<br>
*>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors<br>
*>          of the symmetric tridiagonal matrix.<br>
*>          If  COMPZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1.<br>
*>          If eigenvectors are desired, then LDZ >= max(1,N).<br>
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
*>          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.<br>
*>          If COMPZ = 'V' and N > 1, LWORK must be at least N*N.<br>
*>          Note that for COMPZ = 'V', then if N is less than or<br>
*>          equal to the minimum divide size, usually 25, then LWORK need<br>
*>          only be 1.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK, RWORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array,<br>
*>                                         dimension (LRWORK)<br>
*>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of the array RWORK.<br>
*>          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.<br>
*>          If COMPZ = 'V' and N > 1, LRWORK must be at least<br>
*>                         1 + 3*N + 2*N*lg N + 4*N**2 ,<br>
*>                         where lg( N ) = smallest integer k such<br>
*>                         that 2**k >= N.<br>
*>          If COMPZ = 'I' and N > 1, LRWORK must be at least<br>
*>                         1 + 4*N + 2*N**2 .<br>
*>          Note that for COMPZ = 'I' or 'V', then if N is less than or<br>
*>          equal to the minimum divide size, usually 25, then LRWORK<br>
*>          need only be max(1,2*(N-1)).<br>
*><br>
*>          If LRWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK, RWORK<br>
*>          and IWORK arrays, returns these values as the first entries<br>
*>          of the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
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
*>          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.<br>
*>          If COMPZ = 'V' or N > 1,  LIWORK must be at least<br>
*>                                    6 + 6*N + 5*N*lg N.<br>
*>          If COMPZ = 'I' or N > 1,  LIWORK must be at least<br>
*>                                    3 + 5*N .<br>
*>          Note that for COMPZ = 'I' or 'V', then if N is less than or<br>
*>          equal to the minimum divide size, usually 25, then LIWORK<br>
*>          need only be 1.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK, RWORK<br>
*>          and IWORK arrays, returns these values as the first entries<br>
*>          of the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void zstedc_(CHARACTER COMPZ,INTEGER N,double[] D,double[] E,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZSTEGR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSTEGR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstegr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstegr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstegr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,<br>
*                  ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,<br>
*                  LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE<br>
*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N<br>
*       DOUBLE PRECISION ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )<br>
*       COMPLEX*16         Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSTEGR computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has<br>
*> a well defined set of pairwise different real eigenvalues, the corresponding<br>
*> real eigenvectors are pairwise orthogonal.<br>
*><br>
*> The spectrum may be computed either completely or partially by specifying<br>
*> either an interval (VL,VU] or a range of indices IL:IU for the desired<br>
*> eigenvalues.<br>
*><br>
*> ZSTEGR is a compatibility wrapper around the improved ZSTEMR routine.<br>
*> See DSTEMR for further details.<br>
*><br>
*> One important change is that the ABSTOL parameter no longer provides any<br>
*> benefit and hence is no longer used.<br>
*><br>
*> Note : ZSTEGR and ZSTEMR work only on machines which follow<br>
*> IEEE-754 floating-point standard in their handling of infinities and<br>
*> NaNs.  Normal execution may create these exceptiona values and hence<br>
*> may abort due to a floating point exception in environments which<br>
*> do not conform to the IEEE-754 standard.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': all eigenvalues will be found.<br>
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]<br>
*>                 will be found.<br>
*>          = 'I': the IL-th through IU-th eigenvalues will be found.<br>
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
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the N diagonal elements of the tridiagonal matrix<br>
*>          T. On exit, D is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the (N-1) subdiagonal elements of the tridiagonal<br>
*>          matrix T in elements 1 to N-1 of E. E(N) need not be set on<br>
*>          input, but is used internally as workspace.<br>
*>          On exit, E is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
*><br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          Unused.  Was the absolute error tolerance for the<br>
*>          eigenvalues/eigenvectors in previous versions.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The total number of eigenvalues found.  0 <= M <= N.<br>
*>          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M) )<br>
*>          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix T<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*>          Supplying N columns is always safe.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', then LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISUPPZ<br>
*> \verbatim<br>
*>          ISUPPZ is INTEGER ARRAY, dimension ( 2*max(1,M) )<br>
*>          The support of the eigenvectors in Z, i.e., the indices<br>
*>          indicating the nonzero elements in Z. The i-th computed eigenvector<br>
*>          is nonzero only in elements ISUPPZ( 2*i-1 ) through<br>
*>          ISUPPZ( 2*i ). This is relevant in the case when the matrix<br>
*>          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal<br>
*>          (and minimal) LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,18*N)<br>
*>          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.<br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (LIWORK)<br>
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.  LIWORK >= max(1,10*N)<br>
*>          if the eigenvectors are desired, and LIWORK >= max(1,8*N)<br>
*>          if only the eigenvalues are to be computed.<br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal size of the IWORK array,<br>
*>          returns this value as the first entry of the IWORK array, and<br>
*>          no error message related to LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          On exit, INFO<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = 1X, internal error in DLARRE,<br>
*>                if INFO = 2X, internal error in ZLARRV.<br>
*>                Here, the digit X = ABS( IINFO ) < 10, where IINFO is<br>
*>                the nonzero error code returned by DLARRE or<br>
*>                ZLARRV, respectively.<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Inderjit Dhillon, IBM Almaden, USA \n<br>
*> Osni Marques, LBNL/NERSC, USA \n<br>
*> Christof Voemel, LBNL/NERSC, USA \n<br>
*<br>
*  =====================================================================<br>
*/
	public void zstegr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,double[] D,double[] E,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,int[] ISUPPZ,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZSTEIN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSTEIN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstein.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstein.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstein.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,<br>
*                          IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDZ, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),<br>
*      $                   IWORK( * )<br>
*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )<br>
*       COMPLEX*16         Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSTEIN computes the eigenvectors of a real symmetric tridiagonal<br>
*> matrix T corresponding to specified eigenvalues, using inverse<br>
*> iteration.<br>
*><br>
*> The maximum number of iterations allowed for each eigenvector is<br>
*> specified by an internal parameter MAXITS (currently set to 5).<br>
*><br>
*> Although the eigenvectors are real, they are stored in a complex<br>
*> array, which may be passed to ZUNMTR or ZUPMTR for back<br>
*> transformation to the eigenvectors of a complex Hermitian matrix<br>
*> which was reduced to tridiagonal form.<br>
*><br>
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
*> \param[in] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The n diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The (n-1) subdiagonal elements of the tridiagonal matrix<br>
*>          T, stored in elements 1 to N-1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of eigenvectors to be found.  0 <= M <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          The first M elements of W contain the eigenvalues for<br>
*>          which eigenvectors are to be computed.  The eigenvalues<br>
*>          should be grouped by split-off block and ordered from<br>
*>          smallest to largest within the block.  ( The output array<br>
*>          W from DSTEBZ with ORDER = 'B' is expected here. )<br>
*> \endverbatim<br>
*><br>
*> \param[in] IBLOCK<br>
*> \verbatim<br>
*>          IBLOCK is INTEGER array, dimension (N)<br>
*>          The submatrix indices associated with the corresponding<br>
*>          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to<br>
*>          the first submatrix from the top, =2 if W(i) belongs to<br>
*>          the second submatrix, etc.  ( The output array IBLOCK<br>
*>          from DSTEBZ is expected here. )<br>
*> \endverbatim<br>
*><br>
*> \param[in] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into submatrices.<br>
*>          The first submatrix consists of rows/columns 1 to<br>
*>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1<br>
*>          through ISPLIT( 2 ), etc.<br>
*>          ( The output array ISPLIT from DSTEBZ is expected here. )<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, M)<br>
*>          The computed eigenvectors.  The eigenvector associated<br>
*>          with the eigenvalue W(i) is stored in the i-th column of<br>
*>          Z.  Any vector which fails to converge is set to its current<br>
*>          iterate after MAXITS iterations.<br>
*>          The imaginary parts of the eigenvectors are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (M)<br>
*>          On normal exit, all elements of IFAIL are zero.<br>
*>          If one or more eigenvectors fail to converge after<br>
*>          MAXITS iterations, then their indices are stored in<br>
*>          array IFAIL.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, then i eigenvectors failed to converge<br>
*>               in MAXITS iterations.  Their indices are stored in<br>
*>               array IFAIL.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  MAXITS  INTEGER, default = 5<br>
*>          The maximum number of iterations performed.<br>
*><br>
*>  EXTRA   INTEGER, default = 2<br>
*>          The number of iterations performed after norm growth<br>
*>          criterion is satisfied, should be at least 1.<br>
*> \endverbatim<br>
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
	public void zstein_(INTEGER N,double[] D,double[] E,INTEGER M,double[] W,int[] IBLOCK,int[] ISPLIT,double[] Z,INTEGER LDZ,double[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZSTEMR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSTEMR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstemr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstemr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstemr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,<br>
*                          M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK,<br>
*                          IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE<br>
*       LOGICAL            TRYRAC<br>
*       INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N<br>
*       DOUBLE PRECISION VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )<br>
*       COMPLEX*16         Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSTEMR computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has<br>
*> a well defined set of pairwise different real eigenvalues, the corresponding<br>
*> real eigenvectors are pairwise orthogonal.<br>
*><br>
*> The spectrum may be computed either completely or partially by specifying<br>
*> either an interval (VL,VU] or a range of indices IL:IU for the desired<br>
*> eigenvalues.<br>
*><br>
*> Depending on the number of desired eigenvalues, these are computed either<br>
*> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are<br>
*> computed by the use of various suitable L D L^T factorizations near clusters<br>
*> of close eigenvalues (referred to as RRRs, Relatively Robust<br>
*> Representations). An informal sketch of the algorithm follows.<br>
*><br>
*> For each unreduced block (submatrix) of T,<br>
*>    (a) Compute T - sigma I  = L D L^T, so that L and D<br>
*>        define all the wanted eigenvalues to high relative accuracy.<br>
*>        This means that small relative changes in the entries of D and L<br>
*>        cause only small relative changes in the eigenvalues and<br>
*>        eigenvectors. The standard (unfactored) representation of the<br>
*>        tridiagonal matrix T does not have this property in general.<br>
*>    (b) Compute the eigenvalues to suitable accuracy.<br>
*>        If the eigenvectors are desired, the algorithm attains full<br>
*>        accuracy of the computed eigenvalues only right before<br>
*>        the corresponding vectors have to be computed, see steps c) and d).<br>
*>    (c) For each cluster of close eigenvalues, select a new<br>
*>        shift close to the cluster, find a new factorization, and refine<br>
*>        the shifted eigenvalues to suitable accuracy.<br>
*>    (d) For each eigenvalue with a large enough relative separation compute<br>
*>        the corresponding eigenvector by forming a rank revealing twisted<br>
*>        factorization. Go back to (c) for any clusters that remain.<br>
*><br>
*> For more details, see:<br>
*> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations<br>
*>   to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"<br>
*>   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.<br>
*> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and<br>
*>   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,<br>
*>   2004.  Also LAPACK Working Note 154.<br>
*> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric<br>
*>   tridiagonal eigenvalue/eigenvector problem",<br>
*>   Computer Science Division Technical Report No. UCB/CSD-97-971,<br>
*>   UC Berkeley, May 1997.<br>
*><br>
*> Further Details<br>
*> 1.ZSTEMR works only on machines which follow IEEE-754<br>
*> floating-point standard in their handling of infinities and NaNs.<br>
*> This permits the use of efficient inner loops avoiding a check for<br>
*> zero divisors.<br>
*><br>
*> 2. LAPACK routines can be used to reduce a complex Hermitean matrix to<br>
*> real symmetric tridiagonal form.<br>
*><br>
*> (Any complex Hermitean tridiagonal matrix has real values on its diagonal<br>
*> and potentially complex numbers on its off-diagonals. By applying a<br>
*> similarity transform with an appropriate diagonal matrix<br>
*> diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex Hermitean<br>
*> matrix can be transformed into a real symmetric matrix and complex<br>
*> arithmetic can be entirely avoided.)<br>
*><br>
*> While the eigenvectors of the real symmetric tridiagonal matrix are real,<br>
*> the eigenvectors of original complex Hermitean matrix have complex entries<br>
*> in general.<br>
*> Since LAPACK drivers overwrite the matrix data with the eigenvectors,<br>
*> ZSTEMR accepts complex workspace to facilitate interoperability<br>
*> with ZUNMTR or ZUPMTR.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': all eigenvalues will be found.<br>
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]<br>
*>                 will be found.<br>
*>          = 'I': the IL-th through IU-th eigenvalues will be found.<br>
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
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the N diagonal elements of the tridiagonal matrix<br>
*>          T. On exit, D is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the (N-1) subdiagonal elements of the tridiagonal<br>
*>          matrix T in elements 1 to N-1 of E. E(N) need not be set on<br>
*>          input, but is used internally as workspace.<br>
*>          On exit, E is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
*><br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The total number of eigenvalues found.  0 <= M <= N.<br>
*>          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M) )<br>
*>          If JOBZ = 'V', and if INFO = 0, then the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix T<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and can be computed with a workspace<br>
*>          query by setting NZC = -1, see below.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', then LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] NZC<br>
*> \verbatim<br>
*>          NZC is INTEGER<br>
*>          The number of eigenvectors to be held in the array Z.<br>
*>          If RANGE = 'A', then NZC >= max(1,N).<br>
*>          If RANGE = 'V', then NZC >= the number of eigenvalues in (VL,VU].<br>
*>          If RANGE = 'I', then NZC >= IU-IL+1.<br>
*>          If NZC = -1, then a workspace query is assumed; the<br>
*>          routine calculates the number of columns of the array Z that<br>
*>          are needed to hold the eigenvectors.<br>
*>          This value is returned as the first entry of the Z array, and<br>
*>          no error message related to NZC is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISUPPZ<br>
*> \verbatim<br>
*>          ISUPPZ is INTEGER ARRAY, dimension ( 2*max(1,M) )<br>
*>          The support of the eigenvectors in Z, i.e., the indices<br>
*>          indicating the nonzero elements in Z. The i-th computed eigenvector<br>
*>          is nonzero only in elements ISUPPZ( 2*i-1 ) through<br>
*>          ISUPPZ( 2*i ). This is relevant in the case when the matrix<br>
*>          is split. ISUPPZ is only accessed when JOBZ is 'V' and N > 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] TRYRAC<br>
*> \verbatim<br>
*>          TRYRAC is LOGICAL<br>
*>          If TRYRAC.EQ..TRUE., indicates that the code should check whether<br>
*>          the tridiagonal matrix defines its eigenvalues to high relative<br>
*>          accuracy.  If so, the code uses relative-accuracy preserving<br>
*>          algorithms that might be (a bit) slower depending on the matrix.<br>
*>          If the matrix does not define its eigenvalues to high relative<br>
*>          accuracy, the code can uses possibly faster algorithms.<br>
*>          If TRYRAC.EQ..FALSE., the code is not required to guarantee<br>
*>          relatively accurate eigenvalues and can use the fastest possible<br>
*>          techniques.<br>
*>          On exit, a .TRUE. TRYRAC will be set to .FALSE. if the matrix<br>
*>          does not define its eigenvalues to high relative accuracy.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal<br>
*>          (and minimal) LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK. LWORK >= max(1,18*N)<br>
*>          if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.<br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (LIWORK)<br>
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.  LIWORK >= max(1,10*N)<br>
*>          if the eigenvectors are desired, and LIWORK >= max(1,8*N)<br>
*>          if only the eigenvalues are to be computed.<br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal size of the IWORK array,<br>
*>          returns this value as the first entry of the IWORK array, and<br>
*>          no error message related to LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          On exit, INFO<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = 1X, internal error in DLARRE,<br>
*>                if INFO = 2X, internal error in ZLARRV.<br>
*>                Here, the digit X = ABS( IINFO ) < 10, where IINFO is<br>
*>                the nonzero error code returned by DLARRE or<br>
*>                ZLARRV, respectively.<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Beresford Parlett, University of California, Berkeley, USA \n<br>
*> Jim Demmel, University of California, Berkeley, USA \n<br>
*> Inderjit Dhillon, University of Texas, Austin, USA \n<br>
*> Osni Marques, LBNL/NERSC, USA \n<br>
*> Christof Voemel, University of California, Berkeley, USA \n<br>
*<br>
*  =====================================================================<br>
*/
	public void zstemr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,double[] D,double[] E,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,INTEGER M,double[] W,double[] Z,INTEGER LDZ,INTEGER NZC,int[] ISUPPZ,LOGICAL TRYRAC,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZSTEQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSTEQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsteqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsteqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsteqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPZ<br>
*       INTEGER            INFO, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   D( * ), E( * ), WORK( * )<br>
*       COMPLEX*16         Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a<br>
*> symmetric tridiagonal matrix using the implicit QL or QR method.<br>
*> The eigenvectors of a full or band complex Hermitian matrix can also<br>
*> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this<br>
*> matrix to tridiagonal form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only.<br>
*>          = 'V':  Compute eigenvalues and eigenvectors of the original<br>
*>                  Hermitian matrix.  On entry, Z must contain the<br>
*>                  unitary matrix used to reduce the original matrix<br>
*>                  to tridiagonal form.<br>
*>          = 'I':  Compute eigenvalues and eigenvectors of the<br>
*>                  tridiagonal matrix.  Z is initialized to the identity<br>
*>                  matrix.<br>
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
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          On entry, the diagonal elements of the tridiagonal matrix.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix.<br>
*>          On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          On entry, if  COMPZ = 'V', then Z contains the unitary<br>
*>          matrix used in the reduction to tridiagonal form.<br>
*>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the<br>
*>          orthonormal eigenvectors of the original Hermitian matrix,<br>
*>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors<br>
*>          of the symmetric tridiagonal matrix.<br>
*>          If COMPZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          eigenvectors are desired, then  LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2))<br>
*>          If COMPZ = 'N', then WORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  the algorithm has failed to find all the eigenvalues in<br>
*>                a total of 30*N iterations; if INFO = i, then i<br>
*>                elements of E have not converged to zero; on exit, D<br>
*>                and E contain the elements of a symmetric tridiagonal<br>
*>                matrix which is unitarily similar to the original<br>
*>                matrix.<br>
*> \endverbatim<br>
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
	public void zsteqr_(CHARACTER COMPZ,INTEGER N,double[] D,double[] E,double[] Z,INTEGER LDZ,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsycon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsycon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsycon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYCON estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a complex symmetric matrix A using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by ZSYTRF.<br>
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
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSYTRF.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is DOUBLE PRECISION<br>
*>          The 1-norm of the original matrix A.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsycon_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYCONV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYCONV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyconv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyconv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyconv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, WAY<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYCONV converts A given by ZHETRF into L and D or vice-versa.<br>
*> Get nondiagonal elements of D (returned in workspace) and <br>
*> apply or reverse permutation done in TRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] WAY<br>
*> \verbatim<br>
*>          WAY is CHARACTER*1<br>
*>          = 'C': Convert <br>
*>          = 'R': Revert<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSYTRF.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is COMPLEX*16 array, dimension (N)<br>
*>          E stores the supdiagonal/subdiagonal of the symmetric 1-by-1<br>
*>          or 2-by-2 block diagonal matrix D in LDLT.<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsyconv_(CHARACTER UPLO,CHARACTER WAY,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] E,INTEGER INFO);
/**
*> \brief \b ZSYCON_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYCON_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsycon_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsycon_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsycon_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYCON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND,<br>
*                               WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYCON_ROOK estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a complex symmetric matrix A using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by ZSYTRF_ROOK.<br>
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
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSYTRF_ROOK.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF_ROOK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ANORM<br>
*> \verbatim<br>
*>          ANORM is DOUBLE PRECISION<br>
*>          The 1-norm of the original matrix A.<br>
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
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*> \verbatim<br>
*><br>
*>   November 2015, Igor Kozachenko,<br>
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
	public void zsycon_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYEQUB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYEQUB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyequb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyequb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyequb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, N<br>
*       DOUBLE PRECISION   AMAX, SCOND<br>
*       CHARACTER          UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       DOUBLE PRECISION   S( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYEQUB computes row and column scalings intended to equilibrate a<br>
*> symmetric matrix A and reduce its condition number<br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The N-by-N symmetric matrix whose scaling<br>
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
*>          S is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, S contains the scale factors for A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SCOND<br>
*> \verbatim<br>
*>          SCOND is DOUBLE PRECISION<br>
*>          If INFO = 0, S contains the ratio of the smallest S(i) to<br>
*>          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too<br>
*>          large nor too small, it is not worth scaling by S.<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (3*N)<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n<br>
*>  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n<br>
*>  DOI 10.1023/B:NUMA.0000016606.32820.69 \n <br>
*>  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf<br>
*><br>
*  =====================================================================<br>
*/
	public void zsyequb_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYMV computes a matrix-vector product for a complex symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYMV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsymv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsymv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsymv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INCX, INCY, LDA, N<br>
*       COMPLEX*16         ALPHA, BETA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), X( * ), Y( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n symmetric matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the array A is to be referenced as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension ( LDA, N )<br>
*>           Before entry, with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced.<br>
*>           Before entry, with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           upper triangular part of A is not referenced.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, N ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the N-<br>
*>           element vector x.<br>
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
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y. On exit, Y is overwritten by the updated<br>
*>           vector y.<br>
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
*> \ingroup complex16SYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void zsymv_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZSYR performs the symmetric rank-1 update of a complex symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INCX, LDA, N<br>
*       COMPLEX*16         ALPHA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), X( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYR   performs the symmetric rank 1 operation<br>
*><br>
*>    A := alpha*x*x**H + A,<br>
*><br>
*> where alpha is a complex scalar, x is an n element vector and A is an<br>
*> n by n symmetric matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the array A is to be referenced as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of A<br>
*>                                  is to be referenced.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the N-<br>
*>           element vector x.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension ( LDA, N )<br>
*>           Before entry, with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry, with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           upper triangular part of A is not referenced. On exit, the<br>
*>           lower triangular part of the array A is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, N ).<br>
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
*> \ingroup complex16SYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void zsyr_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] A,INTEGER LDA);
/**
*> \brief \b ZSYRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,<br>
*                          X, LDX, FERR, BERR, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is symmetric indefinite, and<br>
*> provides error bounds and backward error estimates for the solution.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N<br>
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
*>          AF is COMPLEX*16 array, dimension (LDAF,N)<br>
*>          The factored form of the matrix A.  AF contains the block<br>
*>          diagonal matrix D and the multipliers used to obtain the<br>
*>          factor U or L from the factorization A = U*D*U**T or<br>
*>          A = L*D*L**T as computed by ZSYTRF.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF.<br>
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
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by ZSYTRS.<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsyrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZSYRFSX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYRFSX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyrfsx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyrfsx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyrfsx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYRFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                           S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS,<br>
*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,<br>
*                           WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, EQUED<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX, * ), WORK( * )<br>
*       DOUBLE PRECISION   S( * ), PARAMS( * ), BERR( * ), RWORK( * ),<br>
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
*>    ZSYRFSX improves the computed solution to a system of linear<br>
*>    equations when the coefficient matrix is symmetric indefinite, and<br>
*>    provides error bounds and backward error estimates for the<br>
*>    solution.  In addition to normwise error bound, the code provides<br>
*>    maximum componentwise error bound if possible.  See comments for<br>
*>    ERR_BNDS_NORM and ERR_BNDS_COMP for details of the error bounds.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>     The symmetric matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>     upper triangular part of A contains the upper triangular<br>
*>     part of the matrix A, and the strictly lower triangular<br>
*>     part of A is not referenced.  If UPLO = 'L', the leading<br>
*>     N-by-N lower triangular part of A contains the lower<br>
*>     triangular part of the matrix A, and the strictly upper<br>
*>     triangular part of A is not referenced.<br>
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
*>          AF is COMPLEX*16 array, dimension (LDAF,N)<br>
*>     The factored form of the matrix A.  AF contains the block<br>
*>     diagonal matrix D and the multipliers used to obtain the<br>
*>     factor U or L from the factorization A = U*D*U**T or A =<br>
*>     L*D*L**T as computed by DSYTRF.<br>
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
*>     as determined by DSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION array, dimension (N)<br>
*>     The scale factors for A.  If EQUED = 'Y', A is multiplied on<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
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
*>          PARAMS is DOUBLE PRECISION array, dimension NPARAMS<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsyrfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZSYSV computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,<br>
*                         LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS<br>
*> matrices.<br>
*><br>
*> The diagonal pivoting method is used to factor A as<br>
*>    A = U * D * U**T,  if UPLO = 'U', or<br>
*>    A = L * D * L**T,  if UPLO = 'L',<br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is symmetric and block diagonal with<br>
*> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then<br>
*> used to solve the system of equations A * X = B.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the block diagonal matrix D and the<br>
*>          multipliers used to obtain the factor U or L from the<br>
*>          factorization A = U*D*U**T or A = L*D*L**T as computed by<br>
*>          ZSYTRF.<br>
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
*>          Details of the interchanges and the block structure of D, as<br>
*>          determined by ZSYTRF.  If IPIV(k) > 0, then rows and columns<br>
*>          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1<br>
*>          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,<br>
*>          then rows and columns k-1 and -IPIV(k) were interchanged and<br>
*>          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and<br>
*>          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and<br>
*>          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2<br>
*>          diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >= 1, and for best performance<br>
*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for<br>
*>          ZSYTRF.<br>
*>          for LWORK < N, TRS will be done with Level BLAS 2<br>
*>          for LWORK >= N, TRS will be done with Level BLAS 3<br>
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
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular, so the solution could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zsysv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> ZSYSVX computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,<br>
*                          LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK,<br>
*                          RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          FACT, UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS<br>
*       DOUBLE PRECISION   RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * )<br>
*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYSVX uses the diagonal pivoting factorization to compute the<br>
*> solution to a complex system of linear equations A * X = B,<br>
*> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS<br>
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
*> 1. If FACT = 'N', the diagonal pivoting method is used to factor A.<br>
*>    The form of the factorization is<br>
*>       A = U * D * U**T,  if UPLO = 'U', or<br>
*>       A = L * D * L**T,  if UPLO = 'L',<br>
*>    where U (or L) is a product of permutation and unit upper (lower)<br>
*>    triangular matrices, and D is symmetric and block diagonal with<br>
*>    1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*> 2. If some D(i,i)=0, so that D is exactly singular, then the routine<br>
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
*>          = 'F':  On entry, AF and IPIV contain the factored form<br>
*>                  of A.  A, AF and IPIV will not be modified.<br>
*>          = 'N':  The matrix A will be copied to AF and factored.<br>
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
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N<br>
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
*> \param[in,out] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX*16 array, dimension (LDAF,N)<br>
*>          If FACT = 'F', then AF is an input argument and on entry<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by ZSYTRF.<br>
*><br>
*>          If FACT = 'N', then AF is an output argument and on exit<br>
*>          returns the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T.<br>
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
*>          contains details of the interchanges and the block structure<br>
*>          of D, as determined by ZSYTRF.<br>
*>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>          interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and<br>
*>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =<br>
*>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were<br>
*>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*><br>
*>          If FACT = 'N', then IPIV is an output argument and on exit<br>
*>          contains details of the interchanges and the block structure<br>
*>          of D, as determined by ZSYTRF.<br>
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
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >= max(1,2*N), and for best<br>
*>          performance, when FACT = 'N', LWORK >= max(1,2*N,N*NB), where<br>
*>          NB is the optimal blocksize for ZSYTRF.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, and i is<br>
*>                <= N:  D(i,i) is exactly zero.  The factorization<br>
*>                       has been completed but the factor D is exactly<br>
*>                       singular, so the solution and error bounds could<br>
*>                       not be computed. RCOND = 0 is returned.<br>
*>                = N+1: D is nonsingular, but RCOND is less than machine<br>
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
*> \ingroup complex16SYsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zsysvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZSYSVXX computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYSVXX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysvxx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysvxx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysvxx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                           EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR,<br>
*                           N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP,<br>
*                           NPARAMS, PARAMS, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       DOUBLE PRECISION   RCOND, RPVGRW<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX, * ), WORK( * )<br>
*       DOUBLE PRECISION   S( * ), PARAMS( * ), BERR( * ),<br>
*      $                   ERR_BNDS_NORM( NRHS, * ),<br>
*      $                   ERR_BNDS_COMP( NRHS, * ), RWORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZSYSVXX uses the diagonal pivoting factorization to compute the<br>
*>    solution to a complex*16 system of linear equations A * X = B, where<br>
*>    A is an N-by-N symmetric matrix and X and B are N-by-NRHS<br>
*>    matrices.<br>
*><br>
*>    If requested, both normwise and maximum componentwise error bounds<br>
*>    are returned. ZSYSVXX will return a solution with a tiny<br>
*>    guaranteed error (O(eps) where eps is the working machine<br>
*>    precision) unless the matrix is very ill-conditioned, in which<br>
*>    case a warning is returned. Relevant condition numbers also are<br>
*>    calculated and returned.<br>
*><br>
*>    ZSYSVXX accepts user-provided factorizations and equilibration<br>
*>    factors; see the definitions of the FACT and EQUED options.<br>
*>    Solving with refinement and using a factorization from a previous<br>
*>    ZSYSVXX call will also produce a solution with either O(eps)<br>
*>    errors or warnings, but we cannot make that claim for general<br>
*>    user-provided factorizations and equilibration factors if they<br>
*>    differ from what ZSYSVXX would itself produce.<br>
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
*>      diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B<br>
*><br>
*>    Whether or not the system will be equilibrated depends on the<br>
*>    scaling of the matrix A, but if equilibration is used, A is<br>
*>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.<br>
*><br>
*>    2. If FACT = 'N' or 'E', the LU decomposition is used to factor<br>
*>    the matrix A (after equilibration if FACT = 'E') as<br>
*><br>
*>       A = U * D * U**T,  if UPLO = 'U', or<br>
*>       A = L * D * L**T,  if UPLO = 'L',<br>
*><br>
*>    where U (or L) is a product of permutation and unit upper (lower)<br>
*>    triangular matrices, and D is symmetric and block diagonal with<br>
*>    1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*>    3. If some D(i,i)=0, so that D is exactly singular, then the<br>
*>    routine returns with INFO = i. Otherwise, the factored form of A<br>
*>    is used to estimate the condition number of the matrix A (see<br>
*>    argument RCOND).  If the reciprocal of the condition number is<br>
*>    less than machine precision, the routine still goes on to solve<br>
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
*>    diag(R) so that it solves the original system before<br>
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
*>       = 'F':  On entry, AF and IPIV contain the factored form of A.<br>
*>               If EQUED is not 'N', the matrix A has been<br>
*>               equilibrated with scaling factors given by S.<br>
*>               A, AF, and IPIV are not modified.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>     The symmetric matrix A.  If UPLO = 'U', the leading N-by-N<br>
*>     upper triangular part of A contains the upper triangular<br>
*>     part of the matrix A, and the strictly lower triangular<br>
*>     part of A is not referenced.  If UPLO = 'L', the leading<br>
*>     N-by-N lower triangular part of A contains the lower<br>
*>     triangular part of the matrix A, and the strictly upper<br>
*>     triangular part of A is not referenced.<br>
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
*>          AF is COMPLEX*16 array, dimension (LDAF,N)<br>
*>     If FACT = 'F', then AF is an input argument and on entry<br>
*>     contains the block diagonal matrix D and the multipliers<br>
*>     used to obtain the factor U or L from the factorization A =<br>
*>     U*D*U**T or A = L*D*L**T as computed by DSYTRF.<br>
*><br>
*>     If FACT = 'N', then AF is an output argument and on exit<br>
*>     returns the block diagonal matrix D and the multipliers<br>
*>     used to obtain the factor U or L from the factorization A =<br>
*>     U*D*U**T or A = L*D*L**T.<br>
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
*>     contains details of the interchanges and the block<br>
*>     structure of D, as determined by DSYTRF.  If IPIV(k) > 0,<br>
*>     then rows and columns k and IPIV(k) were interchanged and<br>
*>     D(k,k) is a 1-by-1 diagonal block.  If UPLO = 'U' and<br>
*>     IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and<br>
*>     -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2<br>
*>     diagonal block.  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0,<br>
*>     then rows and columns k+1 and -IPIV(k) were interchanged<br>
*>     and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*><br>
*>     If FACT = 'N', then IPIV is an output argument and on exit<br>
*>     contains details of the interchanges and the block<br>
*>     structure of D, as determined by DSYTRF.<br>
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
*>          S is DOUBLE PRECISION array, dimension (N)<br>
*>     The scale factors for A.  If EQUED = 'Y', A is multiplied on<br>
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
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)<br>
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
*>     for the leading INFO columns of A.<br>
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
*>          PARAMS is DOUBLE PRECISION array, dimension NPARAMS<br>
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
*> \ingroup complex16SYsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zsysvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE ZLA_SYRPVGRW);
/**
*> \brief <b> ZSYSV_ROOK computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYSV_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysv_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysv_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysv_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,<br>
*                         LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYSV_ROOK computes the solution to a complex system of linear<br>
*> equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS<br>
*> matrices.<br>
*><br>
*> The diagonal pivoting method is used to factor A as<br>
*>    A = U * D * U**T,  if UPLO = 'U', or<br>
*>    A = L * D * L**T,  if UPLO = 'L',<br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is symmetric and block diagonal with<br>
*> 1-by-1 and 2-by-2 diagonal blocks.  <br>
*><br>
*> ZSYTRF_ROOK is called to compute the factorization of a complex<br>
*> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal<br>
*> pivoting method.<br>
*><br>
*> The factored form of A is then used to solve the system <br>
*> of equations A * X = B by calling ZSYTRS_ROOK.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the block diagonal matrix D and the<br>
*>          multipliers used to obtain the factor U or L from the<br>
*>          factorization A = U*D*U**T or A = L*D*L**T as computed by<br>
*>          ZSYTRF_ROOK.<br>
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
*>          Details of the interchanges and the block structure of D,<br>
*>          as determined by ZSYTRF_ROOK.<br>
*><br>
*>          If UPLO = 'U':<br>
*>               If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>               were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>               If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>               columns k and -IPIV(k) were interchanged and rows and<br>
*>               columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>               If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>               were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>               If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>               columns k and -IPIV(k) were interchanged and rows and<br>
*>               columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >= 1, and for best performance<br>
*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for<br>
*>          ZSYTRF_ROOK.<br>
*>          <br>
*>          TRS will be done with Level 2 BLAS<br>
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
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular, so the solution could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYsolve<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   November 2015, Igor Kozachenko,<br>
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
	public void zsysv_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZSYSWAPR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYSWAPR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsyswapr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsyswapr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsyswapr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYSWAPR( UPLO, N, A, LDA, I1, I2)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER        UPLO<br>
*       INTEGER          I1, I2, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16       A( LDA, N )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYSWAPR applies an elementary permutation on the rows and the columns of<br>
*> a symmetric matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          On entry, the NB diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by ZSYTRF.<br>
*><br>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original<br>
*>          matrix.  If UPLO = 'U', the upper triangular part of the<br>
*>          inverse is formed and the part of A below the diagonal is not<br>
*>          referenced; if UPLO = 'L' the lower triangular part of the<br>
*>          inverse is formed and the part of A above the diagonal is<br>
*>          not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] I1<br>
*> \verbatim<br>
*>          I1 is INTEGER<br>
*>          Index of the first row to swap<br>
*> \endverbatim<br>
*><br>
*> \param[in] I2<br>
*> \verbatim<br>
*>          I2 is INTEGER<br>
*>          Index of the second row to swap<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void zsyswapr_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER I1,INTEGER I2);
/**
*> \brief \b ZSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal pivoting method (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTF2( UPLO, N, A, LDA, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTF2 computes the factorization of a complex symmetric matrix A<br>
*> using the Bunch-Kaufman diagonal pivoting method:<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, U**T is the transpose of U, and D is symmetric and<br>
*> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L (see below for further details).<br>
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
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) = IPIV(k-1) < 0, then rows and columns<br>
*>             k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>             is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) = IPIV(k+1) < 0, then rows and columns<br>
*>             k+1 and -IPIV(k) were interchanged and D(k:k+1,k:k+1)<br>
*>             is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular, and division by zero will occur if it<br>
*>               is used to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**T, where<br>
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,<br>
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to<br>
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    v    0   )   k-s<br>
*>     U(k) =  (   0    I    0   )   s<br>
*>             (   0    0    I   )   n-k<br>
*>                k-s   s   n-k<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).<br>
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),<br>
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).<br>
*><br>
*>  If UPLO = 'L', then A = L*D*L**T, where<br>
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,<br>
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to<br>
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    0     0   )  k-1<br>
*>     L(k) =  (   0    I     0   )  s<br>
*>             (   0    v     I   )  n-k-s+1<br>
*>                k-1   s  n-k-s+1<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).<br>
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),<br>
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  09-29-06 - patch from<br>
*>    Bobby Cheng, MathWorks<br>
*><br>
*>    Replace l.209 and l.377<br>
*>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN<br>
*>    by<br>
*>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN<br>
*><br>
*>  1-96 - Based on modifications by J. Lewis, Boeing Computer Services<br>
*>         Company<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void zsytf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b ZSYTF2_ROOK computes the factorization of a complex symmetric indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZSYTF2_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytf2_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytf2_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytf2_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTF2_ROOK( UPLO, N, A, LDA, IPIV, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTF2_ROOK computes the factorization of a complex symmetric matrix A<br>
*> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, U**T is the transpose of U, and D is symmetric and<br>
*> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L (see below for further details).<br>
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
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>             were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>             were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -k, the k-th argument had an illegal value<br>
*>          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization<br>
*>               has been completed, but the block diagonal matrix D is<br>
*>               exactly singular, and division by zero will occur if it<br>
*>               is used to solve a system of equations.<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**T, where<br>
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,<br>
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to<br>
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    v    0   )   k-s<br>
*>     U(k) =  (   0    I    0   )   s<br>
*>             (   0    0    I   )   n-k<br>
*>                k-s   s   n-k<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).<br>
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),<br>
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).<br>
*><br>
*>  If UPLO = 'L', then A = L*D*L**T, where<br>
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,<br>
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to<br>
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    0     0   )  k-1<br>
*>     L(k) =  (   0    I     0   )  s<br>
*>             (   0    v     I   )  n-k-s+1<br>
*>                k-1   s  n-k-s+1<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).<br>
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),<br>
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).<br>
*> \endverbatim<br>
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
*>  01-01-96 - Based on modifications by<br>
*>    J. Lewis, Boeing Computer Services Company<br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville abd , USA<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void zsytf2_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b ZSYTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRF computes the factorization of a complex symmetric matrix A<br>
*> using the Bunch-Kaufman diagonal pivoting method.  The form of the<br>
*> factorization is<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is symmetric and block diagonal with<br>
*> with 1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L (see below for further details).<br>
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
*>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>          interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and<br>
*>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)<br>
*>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =<br>
*>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were<br>
*>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
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
*>          The length of WORK.  LWORK >=1.  For best performance<br>
*>          LWORK >= N*NB, where NB is the block size returned by ILAENV.<br>
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
*>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization<br>
*>                has been completed, but the block diagonal matrix D is<br>
*>                exactly singular, and division by zero will occur if it<br>
*>                is used to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**T, where<br>
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,<br>
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to<br>
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    v    0   )   k-s<br>
*>     U(k) =  (   0    I    0   )   s<br>
*>             (   0    0    I   )   n-k<br>
*>                k-s   s   n-k<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).<br>
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),<br>
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).<br>
*><br>
*>  If UPLO = 'L', then A = L*D*L**T, where<br>
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,<br>
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to<br>
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    0     0   )  k-1<br>
*>     L(k) =  (   0    I     0   )  s<br>
*>             (   0    v     I   )  n-k-s+1<br>
*>                k-1   s  n-k-s+1<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).<br>
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),<br>
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zsytrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZSYTRF_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRF_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrf_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrf_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrf_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRF_ROOK computes the factorization of a complex symmetric matrix A<br>
*> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.<br>
*> The form of the factorization is<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is symmetric and block diagonal with<br>
*> 1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*> This is the blocked version of the algorithm, calling Level 3 BLAS.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L (see below for further details).<br>
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
*>               If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>               were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>               If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>               columns k and -IPIV(k) were interchanged and rows and<br>
*>               columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>               If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>               were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>               If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>               columns k and -IPIV(k) were interchanged and rows and<br>
*>               columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)).<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >=1.  For best performance<br>
*>          LWORK >= N*NB, where NB is the block size returned by ILAENV.<br>
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
*>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization<br>
*>                has been completed, but the block diagonal matrix D is<br>
*>                exactly singular, and division by zero will occur if it<br>
*>                is used to solve a system of equations.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**T, where<br>
*>     U = P(n)*U(n)* ... *P(k)U(k)* ...,<br>
*>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to<br>
*>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    v    0   )   k-s<br>
*>     U(k) =  (   0    I    0   )   s<br>
*>             (   0    0    I   )   n-k<br>
*>                k-s   s   n-k<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).<br>
*>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),<br>
*>  and A(k,k), and v overwrites A(1:k-2,k-1:k).<br>
*><br>
*>  If UPLO = 'L', then A = L*D*L**T, where<br>
*>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,<br>
*>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to<br>
*>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1<br>
*>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as<br>
*>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such<br>
*>  that if the diagonal block D(k) is of order s (s = 1 or 2), then<br>
*><br>
*>             (   I    0     0   )  k-1<br>
*>     L(k) =  (   0    I     0   )  s<br>
*>             (   0    v     I   )  n-k-s+1<br>
*>                k-1   s  n-k-s+1<br>
*><br>
*>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).<br>
*>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),<br>
*>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   June 2016, Igor Kozachenko,<br>
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
	public void zsytrf_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZSYTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRI computes the inverse of a complex symmetric indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> ZSYTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          On entry, the block diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by ZSYTRF.<br>
*><br>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original<br>
*>          matrix.  If UPLO = 'U', the upper triangular part of the<br>
*>          inverse is formed and the part of A below the diagonal is not<br>
*>          referenced; if UPLO = 'L' the lower triangular part of the<br>
*>          inverse is formed and the part of A above the diagonal is<br>
*>          not referenced.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its<br>
*>               inverse could not be computed.<br>
*> \endverbatim<br>
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
	public void zsytri_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYTRI2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRI2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRI2 computes the inverse of a COMPLEX*16 symmetric indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> ZSYTRF. ZSYTRI2 sets the LEADING DIMENSION of the workspace<br>
*> before calling ZSYTRI2X that actually computes the inverse.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          On entry, the NB diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by ZSYTRF.<br>
*><br>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original<br>
*>          matrix.  If UPLO = 'U', the upper triangular part of the<br>
*>          inverse is formed and the part of A below the diagonal is not<br>
*>          referenced; if UPLO = 'L' the lower triangular part of the<br>
*>          inverse is formed and the part of A above the diagonal is<br>
*>          not referenced.<br>
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
*>          Details of the interchanges and the NB structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N+NB+1)*(NB+3)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          WORK is size >= (N+NB+1)*(NB+3)<br>
*>          If LDWORK = -1, then a workspace query is assumed; the routine<br>
*>           calculates:<br>
*>              - the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array,<br>
*>              - and no error message related to LDWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its<br>
*>               inverse could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsytri2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZSYTRI2X<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRI2X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri2x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri2x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri2x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( N+NB+1,* )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRI2X computes the inverse of a complex symmetric indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> ZSYTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          On entry, the NNB diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by ZSYTRF.<br>
*><br>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original<br>
*>          matrix.  If UPLO = 'U', the upper triangular part of the<br>
*>          inverse is formed and the part of A below the diagonal is not<br>
*>          referenced; if UPLO = 'L' the lower triangular part of the<br>
*>          inverse is formed and the part of A above the diagonal is<br>
*>          not referenced.<br>
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
*>          Details of the interchanges and the NNB structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N+NNB+1,NNB+3)<br>
*> \endverbatim<br>
*><br>
*> \param[in] NB<br>
*> \verbatim<br>
*>          NB is INTEGER<br>
*>          Block size<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its<br>
*>               inverse could not be computed.<br>
*> \endverbatim<br>
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
	public void zsytri2x_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER NB,INTEGER INFO);
/**
*> \brief \b ZSYTRI_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRI_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytri_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytri_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytri_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRI_ROOK computes the inverse of a complex symmetric<br>
*> matrix A using the factorization A = U*D*U**T or A = L*D*L**T<br>
*> computed by ZSYTRF_ROOK.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          On entry, the block diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by ZSYTRF_ROOK.<br>
*><br>
*>          On exit, if INFO = 0, the (symmetric) inverse of the original<br>
*>          matrix.  If UPLO = 'U', the upper triangular part of the<br>
*>          inverse is formed and the part of A below the diagonal is not<br>
*>          referenced; if UPLO = 'L' the lower triangular part of the<br>
*>          inverse is formed and the part of A above the diagonal is<br>
*>          not referenced.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF_ROOK.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its<br>
*>               inverse could not be computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   November 2015, Igor Kozachenko,<br>
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
	public void zsytri_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRS solves a system of linear equations A*X = B with a complex<br>
*> symmetric matrix A using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by ZSYTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSYTRF.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsytrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZSYTRS2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, <br>
*                           WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16       A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRS2 solves a system of linear equations A*X = B with a real<br>
*> symmetric matrix A using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by ZSYTRF and converted by ZSYCONV.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSYTRF.<br>
*>          Note that A is input / output. This might be counter-intuitive,<br>
*>          and one may think that A is input only. A is input / output. This<br>
*>          is because, at the start of the subroutine, we permute A in a<br>
*>          "better" form and then we permute A back to its original form at<br>
*>          the end.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
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
*> \ingroup complex16SYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zsytrs2_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZSYTRS_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZSYTRS_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsytrs_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsytrs_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsytrs_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYTRS_ROOK solves a system of linear equations A*X = B with<br>
*> a complex symmetric matrix A using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by ZSYTRF_ROOK.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the details of the factorization are stored<br>
*>          as an upper or lower triangular matrix.<br>
*>          = 'U':  Upper triangular, form is A = U*D*U**T;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**T.<br>
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
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by ZSYTRF_ROOK.<br>
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
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZSYTRF_ROOK.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup complex16SYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   November 2015, Igor Kozachenko,<br>
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
	public void zsytrs_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);

}