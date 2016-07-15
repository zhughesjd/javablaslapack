package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.CHARACTER;
import net.joshuahughes.jnablaslapack.REAL;
import net.joshuahughes.jnablaslapack.INTEGER;
import net.joshuahughes.jnablaslapack.DOUBLE;

public interface Blas extends Library
{

	public static Blas instance = (Blas) Native.loadLibrary("libblas",Blas.class);

/**
*> \brief \b CAXPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX CA<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*),CY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CAXPY constant times a vector plus a vector.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void caxpy_(INTEGER N,float[] CA,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY);
/**
*> \brief \b CCOPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*),CY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CCOPY copies a vector x to a vector y.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ccopy_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY);
/**
*> \brief \b CDOTC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*),CY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CDOTC forms the dot product of two complex vectors<br>
*>      CDOTC = X^H * Y<br>
*><br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack,  3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float[] cdotc_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY);
/**
*> \brief \b CDOTU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*),CY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CDOTU forms the dot product of two complex vectors<br>
*>      CDOTU = X^T * Y<br>
*><br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float[] cdotu_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY);
/**
*> \brief \b CGBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER INCX,INCY,KL,KU,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CGBMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or<br>
*><br>
*>    y := alpha*A**H*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n band matrix, with kl sub-diagonals and ku super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>           On entry, KL specifies the number of sub-diagonals of the<br>
*>           matrix A. KL must satisfy  0 .le. KL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>           On entry, KU specifies the number of super-diagonals of the<br>
*>           matrix A. KU must satisfy  0 .le. KU.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading ( kl + ku + 1 ) by n part of the<br>
*>           array A must contain the matrix of coefficients, supplied<br>
*>           column by column, with the leading diagonal of the matrix in<br>
*>           row ( ku + 1 ) of the array, the first super-diagonal<br>
*>           starting at position 2 in row ku, the first sub-diagonal<br>
*>           starting at position 1 in row ( ku + 2 ), and so on.<br>
*>           Elements in the array A that do not correspond to elements<br>
*>           in the band matrix (such as the top left ku by ku triangle)<br>
*>           are not referenced.<br>
*>           The following program segment will transfer a band matrix<br>
*>           from conventional full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    K = KU + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )<br>
*>                       A( K + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( kl + ku + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cgbmv_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b CGEMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,M,N<br>
*       CHARACTER TRANSA,TRANSB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CGEMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*op( A )*op( B ) + beta*C,<br>
*><br>
*> where  op( X ) is one of<br>
*><br>
*>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,<br>
*><br>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )<br>
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n',  op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c',  op( A ) = A**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSB<br>
*> \verbatim<br>
*>          TRANSB is CHARACTER*1<br>
*>           On entry, TRANSB specifies the form of op( B ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSB = 'N' or 'n',  op( B ) = B.<br>
*><br>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.<br>
*><br>
*>              TRANSB = 'C' or 'c',  op( B ) = B**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies  the number  of rows  of the  matrix<br>
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N  specifies the number  of columns of the matrix<br>
*>           op( B ) and the number of columns of the matrix C. N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry,  K  specifies  the number of columns of the matrix<br>
*>           op( A ) and the number of rows of the matrix op( B ). K must<br>
*>           be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.<br>
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by m  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, kb ), where kb is<br>
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.<br>
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  n by k  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then<br>
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at<br>
*>           least  max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n  matrix<br>
*>           ( alpha*op( A )*op( B ) + beta*C ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cgemm_(CHARACTER TRANSA,CHARACTER TRANSB,INTEGER M,INTEGER N,INTEGER K,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CGEMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CGEMV performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or<br>
*><br>
*>    y := alpha*A**H*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cgemv_(CHARACTER TRANS,INTEGER M,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b CGERC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CGERC  performs the rank 1 operation<br>
*><br>
*>    A := alpha*x*y**H + A,<br>
*><br>
*> where alpha is a scalar, x is an m element vector, y is an n element<br>
*> vector and A is an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the m<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients. On exit, A is<br>
*>           overwritten by the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cgerc_(INTEGER M,INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] A,INTEGER LDA);
/**
*> \brief \b CGERU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CGERU  performs the rank 1 operation<br>
*><br>
*>    A := alpha*x*y**T + A,<br>
*><br>
*> where alpha is a scalar, x is an m element vector, y is an n element<br>
*> vector and A is an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the m<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients. On exit, A is<br>
*>           overwritten by the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cgeru_(INTEGER M,INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] A,INTEGER LDA);
/**
*> \brief \b CHBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER INCX,INCY,K,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHBMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n hermitian band matrix, with k super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the band matrix A is being supplied as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   The upper triangular part of A is<br>
*>                                  being supplied.<br>
*><br>
*>              UPLO = 'L' or 'l'   The lower triangular part of A is<br>
*>                                  being supplied.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry, K specifies the number of super-diagonals of the<br>
*>           matrix A. K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the hermitian matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer the upper<br>
*>           triangular part of a hermitian band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the hermitian matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer the lower<br>
*>           triangular part of a hermitian band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void chbmv_(CHARACTER UPLO,INTEGER N,INTEGER K,float[] ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b CHEMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER LDA,LDB,LDC,M,N<br>
*       CHARACTER SIDE,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHEMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*A*B + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*B*A + beta*C,<br>
*><br>
*> where alpha and beta are scalars, A is an hermitian matrix and  B and<br>
*> C are m by n matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE  specifies whether  the  hermitian matrix  A<br>
*>           appears on the  left or right  in the  operation as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,<br>
*><br>
*>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of  the  hermitian  matrix   A  is  to  be<br>
*>           referenced as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of the<br>
*>                                  hermitian matrix is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of the<br>
*>                                  hermitian matrix is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies the number of rows of the matrix  C.<br>
*>           M  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix C.<br>
*>           N  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           m  when  SIDE = 'L' or 'l'  and is n  otherwise.<br>
*>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of<br>
*>           the array  A  must contain the  hermitian matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading m by m upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  hermitian matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  m by m  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  hermitian<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of<br>
*>           the array  A  must contain the  hermitian matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading n by n upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  hermitian matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  n by n  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  hermitian<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Note that the imaginary parts  of the diagonal elements need<br>
*>           not be set, they are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, n ).<br>
*>           Before entry, the leading  m by n part of the array  B  must<br>
*>           contain the matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n updated<br>
*>           matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void chemm_(CHARACTER SIDE,CHARACTER UPLO,INTEGER M,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CHEMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHEMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n hermitian matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           lower triangular part of A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           upper triangular part of A is not referenced.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void chemv_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b CHER<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHER   performs the hermitian rank 1 operation<br>
*><br>
*>    A := alpha*x*x**H + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
*> n by n hermitian matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           upper triangular part of A is not referenced. On exit, the<br>
*>           lower triangular part of the array A is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set, they are assumed to be zero, and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cher_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] A,INTEGER LDA);
/**
*> \brief \b CHER2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHER2  performs the hermitian rank 2 operation<br>
*><br>
*>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an n<br>
*> by n hermitian matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           upper triangular part of A is not referenced. On exit, the<br>
*>           lower triangular part of the array A is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set, they are assumed to be zero, and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cher2_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] A,INTEGER LDA);
/**
*> \brief \b CHER2K<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       REAL BETA<br>
*       INTEGER K,LDA,LDB,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHER2K  performs one of the hermitian rank 2k operations<br>
*><br>
*>    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars with  beta  real,  C is an  n by n<br>
*> hermitian matrix and  A and B  are  n by k matrices in the first case<br>
*> and  k by n  matrices in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'    C := alpha*A*B**H          +<br>
*>                                         conjg( alpha )*B*A**H +<br>
*>                                         beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'    C := alpha*A**H*B          +<br>
*>                                         conjg( alpha )*B**H*A +<br>
*>                                         beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns  of the  matrices  A and B,  and on  entry  with<br>
*>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the<br>
*>           matrices  A and B.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, kb ), where kb is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  k by n  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set,  they are assumed to be zero,  and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*><br>
*>  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.<br>
*>     Ed Anderson, Cray Research Inc.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cher2k_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CHERK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER K,LDA,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHERK  performs one of the hermitian rank k operations<br>
*><br>
*>    C := alpha*A*A**H + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**H*A + beta*C,<br>
*><br>
*> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian<br>
*> matrix and  A  is an  n by k  matrix in the  first case and a  k by n<br>
*> matrix in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the<br>
*>           matrix A.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set,  they are assumed to be zero,  and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*><br>
*>  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.<br>
*>     Ed Anderson, Cray Research Inc.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cherk_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,REAL ALPHA,float[] A,INTEGER LDA,REAL BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CHPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHPMV  performs the matrix-vector operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n hermitian matrix, supplied in packed form.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void chpmv_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] AP,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b CHPR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHPR(UPLO,N,ALPHA,X,INCX,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHPR    performs the hermitian rank 1 operation<br>
*><br>
*>    A := alpha*x*x**H + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
*> n by n hermitian matrix, supplied in packed form.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the hermitian matrix<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void chpr_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] AP);
/**
*> \brief \b CHPR2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CHPR2  performs the hermitian rank 2 operation<br>
*><br>
*>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an<br>
*> n by n hermitian matrix, supplied in packed form.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the hermitian matrix<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void chpr2_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] AP);
/**
*> \brief \b CROTG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CROTG(CA,CB,C,S)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX CA,CB,S<br>
*       REAL C<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CROTG determines a complex Givens rotation.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void crotg_(float[] CA,float[] CB,REAL C,float[] S);
/**
*> \brief \b CSCAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSCAL(N,CA,CX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX CA<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CSCAL scales a vector by a constant.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack,  3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cscal_(INTEGER N,float[] CA,float[] CX,INTEGER INCX);
/**
*> \brief \b CSROT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSROT( N, CX, INCX, CY, INCY, C, S )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER           INCX, INCY, N<br>
*       REAL              C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX           CX( * ), CY( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CSROT applies a plane rotation, where the cos and sin (c and s) are real<br>
*> and the vectors cx and cy are complex.<br>
*> jack dongarra, linpack, 3/11/78.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the vectors cx and cy.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CX<br>
*> \verbatim<br>
*>          CX is COMPLEX array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array CX must contain the n<br>
*>           element vector cx. On exit, CX is overwritten by the updated<br>
*>           vector cx.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           CX. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CY<br>
*> \verbatim<br>
*>          CY is COMPLEX array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array CY must contain the n<br>
*>           element vector cy. On exit, CY is overwritten by the updated<br>
*>           vector cy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           CY. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is REAL<br>
*>           On entry, C specifies the cosine, cos.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is REAL<br>
*>           On entry, S specifies the sine, sin.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void csrot_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY,REAL C,REAL S);
/**
*> \brief \b CSSCAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSSCAL(N,SA,CX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL SA<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CSSCAL scales a complex vector by a real constant.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void csscal_(INTEGER N,REAL SA,float[] CX,INTEGER INCX);
/**
*> \brief \b CSWAP<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*),CY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>   CSWAP interchanges two vectors.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void cswap_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY);
/**
*> \brief \b CSYMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER LDA,LDB,LDC,M,N<br>
*       CHARACTER SIDE,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CSYMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*A*B + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*B*A + beta*C,<br>
*><br>
*> where  alpha and beta are scalars, A is a symmetric matrix and  B and<br>
*> C are m by n matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE  specifies whether  the  symmetric matrix  A<br>
*>           appears on the  left or right  in the  operation as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,<br>
*><br>
*>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of  the  symmetric  matrix   A  is  to  be<br>
*>           referenced as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies the number of rows of the matrix  C.<br>
*>           M  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix C.<br>
*>           N  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           m  when  SIDE = 'L' or 'l'  and is n  otherwise.<br>
*>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading m by m upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  m by m  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading n by n upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  n by n  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, n ).<br>
*>           Before entry, the leading  m by n part of the array  B  must<br>
*>           contain the matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n updated<br>
*>           matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void csymm_(CHARACTER SIDE,CHARACTER UPLO,INTEGER M,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CSYR2K<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CSYR2K  performs one of the symmetric rank 2k operations<br>
*><br>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix<br>
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n<br>
*> matrices in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'    C := alpha*A*B**T + alpha*B*A**T +<br>
*>                                         beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'    C := alpha*A**T*B + alpha*B**T*A +<br>
*>                                         beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns  of the  matrices  A and B,  and on  entry  with<br>
*>           TRANS = 'T' or 't',  K  specifies  the number of rows of the<br>
*>           matrices  A and B.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, kb ), where kb is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  k by n  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void csyr2k_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CSYRK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA,BETA<br>
*       INTEGER K,LDA,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CSYRK  performs one of the symmetric rank k operations<br>
*><br>
*>    C := alpha*A*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix<br>
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix<br>
*> in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'T' or 't',  K  specifies  the number of rows of the<br>
*>           matrix A.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void csyrk_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,float[] ALPHA,float[] A,INTEGER LDA,float[] BETA,float[] C,INTEGER LDC);
/**
*> \brief \b CTBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTBMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular band matrix, with ( k + 1 ) diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**H*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctbmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b CTBSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTBSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )<br>
*> diagonals.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**H*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctbsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b CTPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTPMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix, supplied in packed form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**H*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctpmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] AP,float[] X,INTEGER INCX);
/**
*> \brief \b CTPSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTPSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix, supplied in packed form.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**H*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctpsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] AP,float[] X,INTEGER INCX);
/**
*> \brief \b CTRMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTRMM  performs one of the matrix-matrix operations<br>
*><br>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )<br>
*><br>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from<br>
*>           the left or right as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, k ), where k is m<br>
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, n ).<br>
*>           Before entry,  the leading  m by n part of the array  B must<br>
*>           contain the matrix  B,  and  on exit  is overwritten  by the<br>
*>           transformed matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctrmm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b CTRMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTRMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**H*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctrmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b CTRSM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTRSM  solves one of the matrix equations<br>
*><br>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,<br>
*><br>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.<br>
*><br>
*> The matrix X is overwritten on B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry, SIDE specifies whether op( A ) appears on the left<br>
*>           or right of X as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, k ),<br>
*>           where k is m when SIDE = 'L' or 'l'  <br>
*>             and k is n when SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX array of DIMENSION ( LDB, n ).<br>
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
*> \endverbatim<br>
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
*> \ingroup complex_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctrsm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b CTRSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CTRSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**H*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ctrsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b DASUM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DASUM takes the sum of the absolute values.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double dasum_(INTEGER N,double[] DX,INTEGER INCX);
/**
*> \brief \b DAXPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION DA<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*),DY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DAXPY constant times a vector plus a vector.<br>
*>    uses unrolled loops for increments equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void daxpy_(INTEGER N,DOUBLE DA,double[] DX,INTEGER INCX,double[] DY,INTEGER INCY);
/**
*> \brief \b DCABS1<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DCABS1(Z)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 Z<br>
*       ..<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DCABS1 computes |Re(.)| + |Im(.)| of a double complex number <br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public double dcabs1_(double[] Z);
/**
*> \brief \b DCOPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*),DY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DCOPY copies a vector, x, to a vector, y.<br>
*>    uses unrolled loops for increments equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dcopy_(INTEGER N,double[] DX,INTEGER INCX,double[] DY,INTEGER INCY);
/**
*> \brief \b DDOT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*),DY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DDOT forms the dot product of two vectors.<br>
*>    uses unrolled loops for increments equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double ddot_(INTEGER N,double[] DX,INTEGER INCX,double[] DY,INTEGER INCY);
/**
*> \brief \b DGBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER INCX,INCY,KL,KU,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGBMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n band matrix, with kl sub-diagonals and ku super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>           On entry, KL specifies the number of sub-diagonals of the<br>
*>           matrix A. KL must satisfy  0 .le. KL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>           On entry, KU specifies the number of super-diagonals of the<br>
*>           matrix A. KU must satisfy  0 .le. KU.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading ( kl + ku + 1 ) by n part of the<br>
*>           array A must contain the matrix of coefficients, supplied<br>
*>           column by column, with the leading diagonal of the matrix in<br>
*>           row ( ku + 1 ) of the array, the first super-diagonal<br>
*>           starting at position 2 in row ku, the first sub-diagonal<br>
*>           starting at position 1 in row ( ku + 2 ), and so on.<br>
*>           Elements in the array A that do not correspond to elements<br>
*>           in the band matrix (such as the top left ku by ku triangle)<br>
*>           are not referenced.<br>
*>           The following program segment will transfer a band matrix<br>
*>           from conventional full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    K = KU + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )<br>
*>                       A( K + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( kl + ku + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgbmv_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b DGEMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,M,N<br>
*       CHARACTER TRANSA,TRANSB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*op( A )*op( B ) + beta*C,<br>
*><br>
*> where  op( X ) is one of<br>
*><br>
*>    op( X ) = X   or   op( X ) = X**T,<br>
*><br>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )<br>
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n',  op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSB<br>
*> \verbatim<br>
*>          TRANSB is CHARACTER*1<br>
*>           On entry, TRANSB specifies the form of op( B ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSB = 'N' or 'n',  op( B ) = B.<br>
*><br>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.<br>
*><br>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies  the number  of rows  of the  matrix<br>
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N  specifies the number  of columns of the matrix<br>
*>           op( B ) and the number of columns of the matrix C. N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry,  K  specifies  the number of columns of the matrix<br>
*>           op( A ) and the number of rows of the matrix op( B ). K must<br>
*>           be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.<br>
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by m  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is<br>
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.<br>
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  n by k  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then<br>
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at<br>
*>           least  max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n  matrix<br>
*>           ( alpha*op( A )*op( B ) + beta*C ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgemm_(CHARACTER TRANSA,CHARACTER TRANSB,INTEGER M,INTEGER N,INTEGER K,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE BETA,double[] C,INTEGER LDC);
/**
*> \brief \b DGEMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGEMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dgemv_(CHARACTER TRANS,INTEGER M,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b DGER<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DGER   performs the rank 1 operation<br>
*><br>
*>    A := alpha*x*y**T + A,<br>
*><br>
*> where alpha is a scalar, x is an m element vector, y is an n element<br>
*> vector and A is an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the m<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients. On exit, A is<br>
*>           overwritten by the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dger_(INTEGER M,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] A,INTEGER LDA);
/**
*> \brief \b DNRM2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DNRM2 returns the euclidean norm of a vector via the function<br>
*> name, so that<br>
*><br>
*>    DNRM2 := sqrt( x'*x )<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  -- This version written on 25-October-1982.<br>
*>     Modified on 14-October-1993 to inline the call to DLASSQ.<br>
*>     Sven Hammarling, Nag Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double dnrm2_(INTEGER N,double[] X,INTEGER INCX);
/**
*> \brief \b DROT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION C,S<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*),DY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DROT applies a plane rotation.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void drot_(INTEGER N,double[] DX,INTEGER INCX,double[] DY,INTEGER INCY,DOUBLE C,DOUBLE S);
/**
*> \brief \b DROTG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DROTG(DA,DB,C,S)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION C,DA,DB,S<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DROTG construct givens plane rotation.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void drotg_(DOUBLE DA,DOUBLE DB,DOUBLE C,DOUBLE S);
/**
*> \brief \b DROTM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DPARAM(5),DX(*),DY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX<br>
*><br>
*>    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN<br>
*>    (DY**T)<br>
*><br>
*>    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE<br>
*>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.<br>
*>    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..<br>
*><br>
*>    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0<br>
*><br>
*>      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)<br>
*>    H=(          )    (          )    (          )    (          )<br>
*>      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).<br>
*>    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         number of elements in input vector(s)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DX<br>
*> \verbatim<br>
*>          DX is DOUBLE PRECISION array, dimension N<br>
*>         double precision vector with N elements<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>         storage spacing between elements of DX<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DY<br>
*> \verbatim<br>
*>          DY is DOUBLE PRECISION array, dimension N<br>
*>         double precision vector with N elements<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>         storage spacing between elements of DY<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DPARAM<br>
*> \verbatim<br>
*>          DPARAM is DOUBLE PRECISION array, dimension 5<br>
*>     DPARAM(1)=DFLAG<br>
*>     DPARAM(2)=DH11<br>
*>     DPARAM(3)=DH21<br>
*>     DPARAM(4)=DH12<br>
*>     DPARAM(5)=DH22<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void drotm_(INTEGER N,double[] DX,INTEGER INCX,double[] DY,INTEGER INCY,double[] DPARAM);
/**
*> \brief \b DROTMG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION DD1,DD2,DX1,DY1<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DPARAM(5)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS<br>
*>    THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*>    DY2)**T.<br>
*>    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..<br>
*><br>
*>    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0<br>
*><br>
*>      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)<br>
*>    H=(          )    (          )    (          )    (          )<br>
*>      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).<br>
*>    LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22<br>
*>    RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE<br>
*>    VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)<br>
*><br>
*>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE<br>
*>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE<br>
*>    OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in,out] DD1<br>
*> \verbatim<br>
*>          DD1 is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DD2<br>
*> \verbatim<br>
*>          DD2 is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DX1<br>
*> \verbatim<br>
*>          DX1 is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[in] DY1<br>
*> \verbatim<br>
*>          DY1 is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] DPARAM<br>
*> \verbatim<br>
*>          DPARAM is DOUBLE PRECISION array, dimension 5<br>
*>     DPARAM(1)=DFLAG<br>
*>     DPARAM(2)=DH11<br>
*>     DPARAM(3)=DH21<br>
*>     DPARAM(4)=DH12<br>
*>     DPARAM(5)=DH22<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void drotmg_(DOUBLE DD1,DOUBLE DD2,DOUBLE DX1,DOUBLE DY1,double[] DPARAM);
/**
*> \brief \b DSBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER INCX,INCY,K,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSBMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n symmetric band matrix, with k super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the band matrix A is being supplied as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   The upper triangular part of A is<br>
*>                                  being supplied.<br>
*><br>
*>              UPLO = 'L' or 'l'   The lower triangular part of A is<br>
*>                                  being supplied.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry, K specifies the number of super-diagonals of the<br>
*>           matrix A. K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the symmetric matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer the upper<br>
*>           triangular part of a symmetric band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the symmetric matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer the lower<br>
*>           triangular part of a symmetric band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsbmv_(CHARACTER UPLO,INTEGER N,INTEGER K,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b DSCAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSCAL(N,DA,DX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION DA<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DSCAL scales a vector by a constant.<br>
*>    uses unrolled loops for increment equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dscal_(INTEGER N,DOUBLE DA,double[] DX,INTEGER INCX);
/**
*> \brief \b DSDOT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DSDOT(N,SX,INCX,SY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*    AUTHORS<br>
*    =======<br>
*    Lawson, C. L., (JPL), Hanson, R. J., (SNLA), <br>
*    Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Compute the inner product of two vectors with extended<br>
*> precision accumulation and result.<br>
*><br>
*> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY<br>
*> DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),<br>
*> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is<br>
*> defined in a similar way using INCY.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         number of elements in input vector(s)<br>
*> \endverbatim<br>
*><br>
*> \param[in] SX<br>
*> \verbatim<br>
*>          SX is REAL array, dimension(N)<br>
*>         single precision vector with N elements<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          storage spacing between elements of SX<br>
*> \endverbatim<br>
*><br>
*> \param[in] SY<br>
*> \verbatim<br>
*>          SY is REAL array, dimension(N)<br>
*>         single precision vector with N elements<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>         storage spacing between elements of SY<br>
*> \endverbatim<br>
*><br>
*> \result DSDOT<br>
*> \verbatim<br>
*>          DSDOT is DOUBLE PRECISION<br>
*>         DSDOT  double precision dot product (zero if N.LE.0)<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*> \verbatim<br>
*><br>
*>      <br>
*>  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.<br>
*>  Krogh, Basic linear algebra subprograms for Fortran<br>
*>  usage, Algorithm No. 539, Transactions on Mathematical<br>
*>  Software 5, 3 (September 1979), pp. 308-323.<br>
*><br>
*>  REVISION HISTORY  (YYMMDD)<br>
*><br>
*>  791001  DATE WRITTEN<br>
*>  890831  Modified array declarations.  (WRB)<br>
*>  890831  REVISION DATE from Version 3.2<br>
*>  891214  Prologue converted to Version 4.0 format.  (BAB)<br>
*>  920310  Corrected definition of LX in DESCRIPTION.  (WRB)<br>
*>  920501  Reformatted the REFERENCES section.  (WRB)<br>
*>  070118  Reformat to LAPACK style (JL)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double dsdot_(INTEGER N,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY);
/**
*> \brief \b DSPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSPMV  performs the matrix-vector operation<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dspmv_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] AP,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b DSPR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSPR(UPLO,N,ALPHA,X,INCX,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSPR    performs the symmetric rank 1 operation<br>
*><br>
*>    A := alpha*x*x**T + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the lower triangular part of the<br>
*>           updated matrix.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dspr_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] AP);
/**
*> \brief \b DSPR2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSPR2  performs the symmetric rank 2 operation<br>
*><br>
*>    A := alpha*x*y**T + alpha*y*x**T + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the lower triangular part of the<br>
*>           updated matrix.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dspr2_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] AP);
/**
*> \brief \b DSWAP<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*),DY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    interchanges two vectors.<br>
*>    uses unrolled loops for increments equal one.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dswap_(INTEGER N,double[] DX,INTEGER INCX,double[] DY,INTEGER INCY);
/**
*> \brief \b DSYMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER LDA,LDB,LDC,M,N<br>
*       CHARACTER SIDE,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSYMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*A*B + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*B*A + beta*C,<br>
*><br>
*> where alpha and beta are scalars,  A is a symmetric matrix and  B and<br>
*> C are  m by n matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE  specifies whether  the  symmetric matrix  A<br>
*>           appears on the  left or right  in the  operation as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,<br>
*><br>
*>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of  the  symmetric  matrix   A  is  to  be<br>
*>           referenced as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies the number of rows of the matrix  C.<br>
*>           M  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix C.<br>
*>           N  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is<br>
*>           m  when  SIDE = 'L' or 'l'  and is  n otherwise.<br>
*>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading m by m upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  m by m  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading n by n upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  n by n  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least  max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).<br>
*>           Before entry, the leading  m by n part of the array  B  must<br>
*>           contain the matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n updated<br>
*>           matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsymm_(CHARACTER SIDE,CHARACTER UPLO,INTEGER M,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE BETA,double[] C,INTEGER LDC);
/**
*> \brief \b DSYMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSYMV  performs the matrix-vector  operation<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           upper triangular part of A is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsymv_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b DSYR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSYR(UPLO,N,ALPHA,X,INCX,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSYR   performs the symmetric rank 1 operation<br>
*><br>
*>    A := alpha*x*x**T + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
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
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsyr_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] A,INTEGER LDA);
/**
*> \brief \b DSYR2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSYR2  performs the symmetric rank 2 operation<br>
*><br>
*>    A := alpha*x*y**T + alpha*y*x**T + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an n<br>
*> by n symmetric matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
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
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsyr2_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] A,INTEGER LDA);
/**
*> \brief \b DSYR2K<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSYR2K  performs one of the symmetric rank 2k operations<br>
*><br>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix<br>
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n<br>
*> matrices in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*B**T + alpha*B*A**T +<br>
*>                                        beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*B + alpha*B**T*A +<br>
*>                                        beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*B + alpha*B**T*A +<br>
*>                                        beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns  of the  matrices  A and B,  and on  entry  with<br>
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number<br>
*>           of rows of the matrices  A and B.  K must be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  k by n  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsyr2k_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE BETA,double[] C,INTEGER LDC);
/**
*> \brief \b DSYRK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER K,LDA,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DSYRK  performs one of the symmetric rank k operations<br>
*><br>
*>    C := alpha*A*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix<br>
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix<br>
*> in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number<br>
*>           of rows of the matrix  A.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dsyrk_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,DOUBLE ALPHA,double[] A,INTEGER LDA,DOUBLE BETA,double[] C,INTEGER LDC);
/**
*> \brief \b DTBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTBMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular band matrix, with ( k + 1 ) diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**T*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtbmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b DTBSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTBSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )<br>
*> diagonals.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**T*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtbsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b DTPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTPMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix, supplied in packed form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**T*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtpmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] AP,double[] X,INTEGER INCX);
/**
*> \brief \b DTPSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTPSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix, supplied in packed form.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**T*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is DOUBLE PRECISION array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtpsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] AP,double[] X,INTEGER INCX);
/**
*> \brief \b DTRMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTRMM  performs one of the matrix-matrix operations<br>
*><br>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),<br>
*><br>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from<br>
*>           the left or right as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>           A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m<br>
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).<br>
*>           Before entry,  the leading  m by n part of the array  B must<br>
*>           contain the matrix  B,  and  on exit  is overwritten  by the<br>
*>           transformed matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtrmm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB);
/**
*> \brief \b DTRMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTRMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**T*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtrmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b DTRSM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTRSM  solves one of the matrix equations<br>
*><br>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,<br>
*><br>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T.<br>
*><br>
*> The matrix X is overwritten on B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry, SIDE specifies whether op( A ) appears on the left<br>
*>           or right of X as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, k ),<br>
*>           where k is m when SIDE = 'L' or 'l'  <br>
*>             and k is n when SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ).<br>
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
*> \endverbatim<br>
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
*> \ingroup double_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dtrsm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB);
/**
*> \brief \b DTRSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DTRSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**T*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is DOUBLE PRECISION array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void dtrsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b DZASUM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DZASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and<br>
*>    returns a single precision result.<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double dzasum_(INTEGER N,double[] ZX,INTEGER INCX);
/**
*> \brief \b DZNRM2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DZNRM2 returns the euclidean norm of a vector via the function<br>
*> name, so that<br>
*><br>
*>    DZNRM2 := sqrt( x**H*x )<br>
*> \endverbatim<br>
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
*> \ingroup double_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  -- This version written on 25-October-1982.<br>
*>     Modified on 14-October-1993 to inline the call to ZLASSQ.<br>
*>     Sven Hammarling, Nag Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double dznrm2_(INTEGER N,double[] X,INTEGER INCX);
/**
*> \brief \b ICAMAX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ICAMAX(N,CX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public int icamax_(INTEGER N,float[] CX,INTEGER INCX);
/**
*> \brief \b IDAMAX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION IDAMAX(N,DX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION DX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    IDAMAX finds the index of the first element having maximum absolute value.<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public int idamax_(INTEGER N,double[] DX,INTEGER INCX);
/**
*> \brief \b ISAMAX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ISAMAX(N,SX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ISAMAX finds the index of the first element having maximum absolute value.<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public int isamax_(INTEGER N,float[] SX,INTEGER INCX);
/**
*> \brief \b IZAMAX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION IZAMAX(N,ZX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 1/15/85.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public int izamax_(INTEGER N,double[] ZX,INTEGER INCX);
/**
*> \brief \b LSAME<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       LOGICAL FUNCTION LSAME(CA,CB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER CA,CB<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> LSAME returns .TRUE. if CA is the same letter as CB regardless of<br>
*> case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] CA<br>
*> \verbatim<br>
*>          CA is CHARACTER*1<br>
*> \endverbatim<br>
*><br>
*> \param[in] CB<br>
*> \verbatim<br>
*>          CB is CHARACTER*1<br>
*>          CA and CB specify the single characters to be compared.<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*  =====================================================================<br>
*/
	public void lsame_(CHARACTER CA,CHARACTER CB);
/**
*> \brief \b SASUM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SASUM(N,SX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SASUM takes the sum of the absolute values.<br>
*>    uses unrolled loops for increment equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float sasum_(INTEGER N,float[] SX,INTEGER INCX);
/**
*> \brief \b SAXPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL SA<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SAXPY constant times a vector plus a vector.<br>
*>    uses unrolled loops for increments equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void saxpy_(INTEGER N,REAL SA,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY);
/**
*> \brief \b SCABS1<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SCABS1(Z)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX Z<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SCABS1 computes |Re(.)| + |Im(.)| of a complex number<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public float scabs1_(float[] Z);
/**
*> \brief \b SCASUM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SCASUM(N,CX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX CX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SCASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and<br>
*>    returns a single precision result.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float scasum_(INTEGER N,float[] CX,INTEGER INCX);
/**
*> \brief \b SCNRM2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SCNRM2(N,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SCNRM2 returns the euclidean norm of a vector via the function<br>
*> name, so that<br>
*><br>
*>    SCNRM2 := sqrt( x**H*x )<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  -- This version written on 25-October-1982.<br>
*>     Modified on 14-October-1993 to inline the call to CLASSQ.<br>
*>     Sven Hammarling, Nag Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float scnrm2_(INTEGER N,float[] X,INTEGER INCX);
/**
*> \brief \b SCOPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SCOPY copies a vector, x, to a vector, y.<br>
*>    uses unrolled loops for increments equal to 1.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void scopy_(INTEGER N,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY);
/**
*> \brief \b SDOT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SDOT forms the dot product of two vectors.<br>
*>    uses unrolled loops for increments equal to one.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float sdot_(INTEGER N,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY);
/**
*> \brief \b SDSDOT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL SB<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*    PURPOSE<br>
*    =======<br>
*  <br>
*    Compute the inner product of two vectors with extended<br>
*    precision accumulation.<br>
*  <br>
*    Returns S.P. result with dot product accumulated in D.P.<br>
*    SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),<br>
*    where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is<br>
*    defined in a similar way using INCY.<br>
*  <br>
*    AUTHOR<br>
*    ======<br>
*    Lawson, C. L., (JPL), Hanson, R. J., (SNLA),<br>
*    Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)<br>
*  <br>
*    ARGUMENTS <br>
*    =========<br>
*  <br>
*    N      (input) INTEGER<br>
*           number of elements in input vector(s)<br>
*  <br>
*    SB     (input) REAL<br>
*           single precision scalar to be added to inner product<br>
*  <br>
*    SX     (input) REAL array, dimension (N)<br>
*           single precision vector with N elements<br>
*  <br>
*    INCX   (input) INTEGER<br>
*           storage spacing between elements of SX<br>
*  <br>
*    SY     (input) REAL array, dimension (N)<br>
*           single precision vector with N elements<br>
*  <br>
*    INCY   (input) INTEGER<br>
*           storage spacing between elements of SY<br>
*  <br>
*    SDSDOT (output) REAL<br>
*           single precision dot product (SB if N .LE. 0)<br>
*  <br>
*    Further Details<br>
*    ===============<br>
*  <br>
*    REFERENCES<br>
*  <br>
*    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.<br>
*    Krogh, Basic linear algebra subprograms for Fortran<br>
*    usage, Algorithm No. 539, Transactions on Mathematical<br>
*    Software 5, 3 (September 1979), pp. 308-323.<br>
*  <br>
*    REVISION HISTORY  (YYMMDD)<br>
*        <br>
*    791001  DATE WRITTEN<br>
*    890531  Changed all specific intrinsics to generic.  (WRB)<br>
*    890831  Modified array declarations.  (WRB)<br>
*    890831  REVISION DATE from Version 3.2<br>
*    891214  Prologue converted to Version 4.0 format.  (BAB)<br>
*    920310  Corrected definition of LX in DESCRIPTION.  (WRB)<br>
*    920501  Reformatted the REFERENCES section.  (WRB)<br>
*    070118  Reformat to LAPACK coding style<br>
*  <br>
*    =====================================================================<br>
*  <br>
*       .. Local Scalars ..<br>
*       DOUBLE PRECISION DSDOT<br>
*       INTEGER I,KX,KY,NS<br>
*       ..<br>
*       .. Intrinsic Functions ..<br>
*       INTRINSIC DBLE<br>
*       ..<br>
*       DSDOT = SB<br>
*       IF (N.LE.0) THEN<br>
*          SDSDOT = DSDOT<br>
*          RETURN<br>
*       END IF   <br>
*       IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN<br>
*  <br>
*       Code for equal and positive increments.<br>
*  <br>
*          NS = N*INCX<br>
*          DO I = 1,NS,INCX<br>
*             DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))<br>
*          END DO<br>
*       ELSE<br>
*  <br>
*       Code for unequal or nonpositive increments.<br>
*  <br>
*          KX = 1<br>
*          KY = 1<br>
*          IF (INCX.LT.0) KX = 1 + (1-N)*INCX<br>
*          IF (INCY.LT.0) KY = 1 + (1-N)*INCY<br>
*          DO I = 1,N<br>
*             DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))<br>
*             KX = KX + INCX<br>
*             KY = KY + INCY<br>
*          END DO<br>
*       END IF<br>
*       SDSDOT = DSDOT<br>
*       RETURN<br>
*       END<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public float sdsdot_(INTEGER N,REAL SB,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY);
/**
*> \brief \b SGBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER INCX,INCY,KL,KU,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SGBMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n band matrix, with kl sub-diagonals and ku super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>           On entry, KL specifies the number of sub-diagonals of the<br>
*>           matrix A. KL must satisfy  0 .le. KL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>           On entry, KU specifies the number of super-diagonals of the<br>
*>           matrix A. KU must satisfy  0 .le. KU.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading ( kl + ku + 1 ) by n part of the<br>
*>           array A must contain the matrix of coefficients, supplied<br>
*>           column by column, with the leading diagonal of the matrix in<br>
*>           row ( ku + 1 ) of the array, the first super-diagonal<br>
*>           starting at position 2 in row ku, the first sub-diagonal<br>
*>           starting at position 1 in row ( ku + 2 ), and so on.<br>
*>           Elements in the array A that do not correspond to elements<br>
*>           in the band matrix (such as the top left ku by ku triangle)<br>
*>           are not referenced.<br>
*>           The following program segment will transfer a band matrix<br>
*>           from conventional full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    K = KU + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )<br>
*>                       A( K + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( kl + ku + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sgbmv_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b SGEMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,M,N<br>
*       CHARACTER TRANSA,TRANSB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SGEMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*op( A )*op( B ) + beta*C,<br>
*><br>
*> where  op( X ) is one of<br>
*><br>
*>    op( X ) = X   or   op( X ) = X**T,<br>
*><br>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )<br>
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n',  op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c',  op( A ) = A**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSB<br>
*> \verbatim<br>
*>          TRANSB is CHARACTER*1<br>
*>           On entry, TRANSB specifies the form of op( B ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSB = 'N' or 'n',  op( B ) = B.<br>
*><br>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.<br>
*><br>
*>              TRANSB = 'C' or 'c',  op( B ) = B**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies  the number  of rows  of the  matrix<br>
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N  specifies the number  of columns of the matrix<br>
*>           op( B ) and the number of columns of the matrix C. N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry,  K  specifies  the number of columns of the matrix<br>
*>           op( A ) and the number of rows of the matrix op( B ). K must<br>
*>           be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.<br>
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by m  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array of DIMENSION ( LDB, kb ), where kb is<br>
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.<br>
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  n by k  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then<br>
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at<br>
*>           least  max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n  matrix<br>
*>           ( alpha*op( A )*op( B ) + beta*C ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sgemm_(CHARACTER TRANSA,CHARACTER TRANSB,INTEGER M,INTEGER N,INTEGER K,REAL ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL BETA,float[] C,INTEGER LDC);
/**
*> \brief \b SGEMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SGEMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sgemv_(CHARACTER TRANS,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b SGER<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SGER   performs the rank 1 operation<br>
*><br>
*>    A := alpha*x*y**T + A,<br>
*><br>
*> where alpha is a scalar, x is an m element vector, y is an n element<br>
*> vector and A is an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the m<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients. On exit, A is<br>
*>           overwritten by the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sger_(INTEGER M,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] A,INTEGER LDA);
/**
*> \brief \b SNRM2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       REAL FUNCTION SNRM2(N,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SNRM2 returns the euclidean norm of a vector via the function<br>
*> name, so that<br>
*><br>
*>    SNRM2 := sqrt( x'*x ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  -- This version written on 25-October-1982.<br>
*>     Modified on 14-October-1993 to inline the call to SLASSQ.<br>
*>     Sven Hammarling, Nag Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public float snrm2_(INTEGER N,float[] X,INTEGER INCX);
/**
*> \brief \b SROT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL C,S<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    applies a plane rotation.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void srot_(INTEGER N,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY,REAL C,REAL S);
/**
*> \brief \b SROTG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SROTG(SA,SB,C,S)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL C,S,SA,SB<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    SROTG construct givens plane rotation.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void srotg_(REAL SA,REAL SB,REAL C,REAL S);
/**
*> \brief \b SROTM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SPARAM(5),SX(*),SY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX<br>
*><br>
*>    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN<br>
*>    (SX**T)<br>
*><br>
*>    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE<br>
*>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.<br>
*>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..<br>
*><br>
*>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0<br>
*><br>
*>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)<br>
*>    H=(          )    (          )    (          )    (          )<br>
*>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).<br>
*>    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>         number of elements in input vector(s)<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SX<br>
*> \verbatim<br>
*>          SX is REAL array, dimension N<br>
*>         double precision vector with N elements<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>         storage spacing between elements of SX<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SY<br>
*> \verbatim<br>
*>          SY is REAL array, dimension N<br>
*>         double precision vector with N elements<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>         storage spacing between elements of SY<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SPARAM<br>
*> \verbatim<br>
*>          SPARAM is REAL array, dimension 5<br>
*>     SPARAM(1)=SFLAG<br>
*>     SPARAM(2)=SH11<br>
*>     SPARAM(3)=SH21<br>
*>     SPARAM(4)=SH12<br>
*>     SPARAM(5)=SH22<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void srotm_(INTEGER N,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY,float[] SPARAM);
/**
*> \brief \b SROTMG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL SD1,SD2,SX1,SY1<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SPARAM(5)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS<br>
*>    THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*>    SY2)**T.<br>
*>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..<br>
*><br>
*>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0<br>
*><br>
*>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)<br>
*>    H=(          )    (          )    (          )    (          )<br>
*>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).<br>
*>    LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22<br>
*>    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE<br>
*>    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)<br>
*><br>
*>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE<br>
*>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE<br>
*>    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.<br>
*><br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in,out] SD1<br>
*> \verbatim<br>
*>          SD1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SD2<br>
*> \verbatim<br>
*>          SD2 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SX1<br>
*> \verbatim<br>
*>          SX1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] SY1<br>
*> \verbatim<br>
*>          SY1 is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SPARAM<br>
*> \verbatim<br>
*>          SPARAM is REAL array, dimension 5<br>
*>     SPARAM(1)=SFLAG<br>
*>     SPARAM(2)=SH11<br>
*>     SPARAM(3)=SH21<br>
*>     SPARAM(4)=SH12<br>
*>     SPARAM(5)=SH22<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void srotmg_(REAL SD1,REAL SD2,REAL SX1,REAL SY1,float[] SPARAM);
/**
*> \brief \b SSBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER INCX,INCY,K,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n symmetric band matrix, with k super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the band matrix A is being supplied as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   The upper triangular part of A is<br>
*>                                  being supplied.<br>
*><br>
*>              UPLO = 'L' or 'l'   The lower triangular part of A is<br>
*>                                  being supplied.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry, K specifies the number of super-diagonals of the<br>
*>           matrix A. K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the symmetric matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer the upper<br>
*>           triangular part of a symmetric band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the symmetric matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer the lower<br>
*>           triangular part of a symmetric band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssbmv_(CHARACTER UPLO,INTEGER N,INTEGER K,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b SSCAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSCAL(N,SA,SX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL SA<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    scales a vector by a constant.<br>
*>    uses unrolled loops for increment equal to 1.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sscal_(INTEGER N,REAL SA,float[] SX,INTEGER INCX);
/**
*> \brief \b SSPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPMV  performs the matrix-vector operation<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is REAL array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sspmv_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] AP,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b SSPR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPR    performs the symmetric rank 1 operation<br>
*><br>
*>    A := alpha*x*x**T + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is REAL array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the lower triangular part of the<br>
*>           updated matrix.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sspr_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] AP);
/**
*> \brief \b SSPR2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPR2  performs the symmetric rank 2 operation<br>
*><br>
*>    A := alpha*x*y**T + alpha*y*x**T + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is REAL array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the symmetric matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the lower triangular part of the<br>
*>           updated matrix.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sspr2_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] AP);
/**
*> \brief \b SSWAP<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL SX(*),SY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    interchanges two vectors.<br>
*>    uses unrolled loops for increments equal to 1.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void sswap_(INTEGER N,float[] SX,INTEGER INCX,float[] SY,INTEGER INCY);
/**
*> \brief \b SSYMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER LDA,LDB,LDC,M,N<br>
*       CHARACTER SIDE,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*A*B + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*B*A + beta*C,<br>
*><br>
*> where alpha and beta are scalars,  A is a symmetric matrix and  B and<br>
*> C are  m by n matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE  specifies whether  the  symmetric matrix  A<br>
*>           appears on the  left or right  in the  operation as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,<br>
*><br>
*>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of  the  symmetric  matrix   A  is  to  be<br>
*>           referenced as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies the number of rows of the matrix  C.<br>
*>           M  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix C.<br>
*>           N  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, ka ), where ka is<br>
*>           m  when  SIDE = 'L' or 'l'  and is  n otherwise.<br>
*>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading m by m upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  m by m  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading n by n upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  n by n  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least  max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array of DIMENSION ( LDB, n ).<br>
*>           Before entry, the leading  m by n part of the array  B  must<br>
*>           contain the matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n updated<br>
*>           matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssymm_(CHARACTER SIDE,CHARACTER UPLO,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL BETA,float[] C,INTEGER LDC);
/**
*> \brief \b SSYMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYMV  performs the matrix-vector  operation<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           upper triangular part of A is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssymv_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY);
/**
*> \brief \b SSYR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYR   performs the symmetric rank 1 operation<br>
*><br>
*>    A := alpha*x*x**T + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
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
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssyr_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] A,INTEGER LDA);
/**
*> \brief \b SSYR2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYR2  performs the symmetric rank 2 operation<br>
*><br>
*>    A := alpha*x*y**T + alpha*y*x**T + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an n<br>
*> by n symmetric matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the symmetric matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
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
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssyr2_(CHARACTER UPLO,INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] A,INTEGER LDA);
/**
*> \brief \b SSYR2K<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYR2K  performs one of the symmetric rank 2k operations<br>
*><br>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix<br>
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n<br>
*> matrices in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*B**T + alpha*B*A**T +<br>
*>                                        beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*B + alpha*B**T*A +<br>
*>                                        beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*B + alpha*B**T*A +<br>
*>                                        beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns  of the  matrices  A and B,  and on  entry  with<br>
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number<br>
*>           of rows of the matrices  A and B.  K must be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array of DIMENSION ( LDB, kb ), where kb is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  k by n  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssyr2k_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,REAL ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL BETA,float[] C,INTEGER LDC);
/**
*> \brief \b SSYRK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA,BETA<br>
*       INTEGER K,LDA,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYRK  performs one of the symmetric rank k operations<br>
*><br>
*>    C := alpha*A*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix<br>
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix<br>
*> in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number<br>
*>           of rows of the matrix  A.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssyrk_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,REAL ALPHA,float[] A,INTEGER LDA,REAL BETA,float[] C,INTEGER LDC);
/**
*> \brief \b STBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STBMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular band matrix, with ( k + 1 ) diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**T*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stbmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b STBSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STBSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )<br>
*> diagonals.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**T*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stbsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b STPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix, supplied in packed form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**T*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is REAL array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stpmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] AP,float[] X,INTEGER INCX);
/**
*> \brief \b STPSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STPSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix, supplied in packed form.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**T*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is REAL array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void stpsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] AP,float[] X,INTEGER INCX);
/**
*> \brief \b STRMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRMM  performs one of the matrix-matrix operations<br>
*><br>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),<br>
*><br>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from<br>
*>           the left or right as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, k ), where k is m<br>
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array of DIMENSION ( LDB, n ).<br>
*>           Before entry,  the leading  m by n part of the array  B must<br>
*>           contain the matrix  B,  and  on exit  is overwritten  by the<br>
*>           transformed matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void strmm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b STRMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**T*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void strmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b STRSM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRSM  solves one of the matrix equations<br>
*><br>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,<br>
*><br>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T.<br>
*><br>
*> The matrix X is overwritten on B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry, SIDE specifies whether op( A ) appears on the left<br>
*>           or right of X as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, k ),<br>
*>           where k is m when SIDE = 'L' or 'l'  <br>
*>             and k is n when SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array of DIMENSION ( LDB, n ).<br>
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
*> \endverbatim<br>
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
*> \ingroup single_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void strsm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
/**
*> \brief \b STRSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> STRSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**T*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup single_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void strsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,float[] A,INTEGER LDA,float[] X,INTEGER INCX);
/**
*> \brief \b XERBLA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE XERBLA( SRNAME, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*(*)      SRNAME<br>
*       INTEGER            INFO<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> XERBLA  is an error handler for the LAPACK routines.<br>
*> It is called by an LAPACK routine if an input parameter has an<br>
*> invalid value.  A message is printed and execution stops.<br>
*><br>
*> Installers may consider modifying the STOP statement in order to<br>
*> call system-specific exception-handling facilities.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SRNAME<br>
*> \verbatim<br>
*>          SRNAME is CHARACTER*(*)<br>
*>          The name of the routine which called XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          The position of the invalid parameter in the parameter list<br>
*>          of the calling routine.<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*  =====================================================================<br>
*/
	public void xerbla_(CHARACTER SRNAME,INTEGER INFO);
/**
*> \brief \b XERBLA_ARRAY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE XERBLA_ARRAY(SRNAME_ARRAY, SRNAME_LEN, INFO)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER SRNAME_LEN, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       CHARACTER(1) SRNAME_ARRAY(SRNAME_LEN)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> XERBLA_ARRAY assists other languages in calling XERBLA, the LAPACK<br>
*> and BLAS error handler.  Rather than taking a Fortran string argument<br>
*> as the function's name, XERBLA_ARRAY takes an array of single<br>
*> characters along with the array's length.  XERBLA_ARRAY then copies<br>
*> up to 32 characters of that array into a Fortran string and passes<br>
*> that to XERBLA.  If called with a non-positive SRNAME_LEN,<br>
*> XERBLA_ARRAY will call XERBLA with a string of all blank characters.<br>
*><br>
*> Say some macro or other device makes XERBLA_ARRAY available to C99<br>
*> by a name lapack_xerbla and with a common Fortran calling convention.<br>
*> Then a C99 program could invoke XERBLA via:<br>
*>    {<br>
*>      int flen = strlen(__func__);<br>
*>      lapack_xerbla(__func__, &flen, &info);<br>
*>    }<br>
*><br>
*> Providing XERBLA_ARRAY is not necessary for intercepting LAPACK<br>
*> errors.  XERBLA_ARRAY calls XERBLA.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SRNAME_ARRAY<br>
*> \verbatim<br>
*>          SRNAME_ARRAY is CHARACTER(1) array, dimension (SRNAME_LEN)<br>
*>          The name of the routine which called XERBLA_ARRAY.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SRNAME_LEN<br>
*> \verbatim<br>
*>          SRNAME_LEN is INTEGER<br>
*>          The length of the name in SRNAME_ARRAY.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          The position of the invalid parameter in the parameter list<br>
*>          of the calling routine.<br>
*> \endverbatim<br>
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
*> \ingroup aux_blas<br>
*<br>
*  =====================================================================<br>
*/
	public void xerbla_array_(char[] SRNAME_ARRAY,INTEGER SRNAME_LEN,INTEGER INFO);
/**
*> \brief \b ZAXPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ZA<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*),ZY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZAXPY constant times a vector plus a vector.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zaxpy_(INTEGER N,double[] ZA,double[] ZX,INTEGER INCX,double[] ZY,INTEGER INCY);
/**
*> \brief \b ZCOPY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*),ZY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZCOPY copies a vector, x, to a vector, y.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, linpack, 4/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zcopy_(INTEGER N,double[] ZX,INTEGER INCX,double[] ZY,INTEGER INCY);
/**
*> \brief \b ZDOTC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*),ZY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZDOTC forms the dot product of two complex vectors<br>
*>      ZDOTC = X^H * Y<br>
*><br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double[] zdotc_(INTEGER N,double[] ZX,INTEGER INCX,double[] ZY,INTEGER INCY);
/**
*> \brief \b ZDOTU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*),ZY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZDOTU forms the dot product of two complex vectors<br>
*>      ZDOTU = X^T * Y<br>
*><br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public double[] zdotu_(INTEGER N,double[] ZX,INTEGER INCX,double[] ZY,INTEGER INCY);
/**
*> \brief \b ZDROT<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, INCY, N<br>
*       DOUBLE PRECISION   C, S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         CX( * ), CY( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Applies a plane rotation, where the cos and sin (c and s) are real<br>
*> and the vectors cx and cy are complex.<br>
*> jack dongarra, linpack, 3/11/78.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the vectors cx and cy.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CX<br>
*> \verbatim<br>
*>          CX is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array CX must contain the n<br>
*>           element vector cx. On exit, CX is overwritten by the updated<br>
*>           vector cx.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           CX. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CY<br>
*> \verbatim<br>
*>          CY is COMPLEX*16 array, dimension at least<br>
*>           ( 1 + ( N - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array CY must contain the n<br>
*>           element vector cy. On exit, CY is overwritten by the updated<br>
*>           vector cy.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           CY. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION<br>
*>           On entry, C specifies the cosine, cos.<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is DOUBLE PRECISION<br>
*>           On entry, S specifies the sine, sin.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void zdrot_(INTEGER N,double[] CX,INTEGER INCX,double[] CY,INTEGER INCY,DOUBLE C,DOUBLE S);
/**
*> \brief \b ZDSCAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZDSCAL(N,DA,ZX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION DA<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZDSCAL scales a vector by a constant.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zdscal_(INTEGER N,DOUBLE DA,double[] ZX,INTEGER INCX);
/**
*> \brief \b ZGBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER INCX,INCY,KL,KU,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZGBMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or<br>
*><br>
*>    y := alpha*A**H*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n band matrix, with kl sub-diagonals and ku super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KL<br>
*> \verbatim<br>
*>          KL is INTEGER<br>
*>           On entry, KL specifies the number of sub-diagonals of the<br>
*>           matrix A. KL must satisfy  0 .le. KL.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KU<br>
*> \verbatim<br>
*>          KU is INTEGER<br>
*>           On entry, KU specifies the number of super-diagonals of the<br>
*>           matrix A. KU must satisfy  0 .le. KU.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading ( kl + ku + 1 ) by n part of the<br>
*>           array A must contain the matrix of coefficients, supplied<br>
*>           column by column, with the leading diagonal of the matrix in<br>
*>           row ( ku + 1 ) of the array, the first super-diagonal<br>
*>           starting at position 2 in row ku, the first sub-diagonal<br>
*>           starting at position 1 in row ( ku + 2 ), and so on.<br>
*>           Elements in the array A that do not correspond to elements<br>
*>           in the band matrix (such as the top left ku by ku triangle)<br>
*>           are not referenced.<br>
*>           The following program segment will transfer a band matrix<br>
*>           from conventional full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    K = KU + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )<br>
*>                       A( K + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( kl + ku + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zgbmv_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZGEMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,M,N<br>
*       CHARACTER TRANSA,TRANSB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZGEMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*op( A )*op( B ) + beta*C,<br>
*><br>
*> where  op( X ) is one of<br>
*><br>
*>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,<br>
*><br>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )<br>
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n',  op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c',  op( A ) = A**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSB<br>
*> \verbatim<br>
*>          TRANSB is CHARACTER*1<br>
*>           On entry, TRANSB specifies the form of op( B ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSB = 'N' or 'n',  op( B ) = B.<br>
*><br>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.<br>
*><br>
*>              TRANSB = 'C' or 'c',  op( B ) = B**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies  the number  of rows  of the  matrix<br>
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N  specifies the number  of columns of the matrix<br>
*>           op( B ) and the number of columns of the matrix C. N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry,  K  specifies  the number of columns of the matrix<br>
*>           op( A ) and the number of rows of the matrix op( B ). K must<br>
*>           be at least  zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.<br>
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by m  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is<br>
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.<br>
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  n by k  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then<br>
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at<br>
*>           least  max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n  matrix<br>
*>           ( alpha*op( A )*op( B ) + beta*C ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zgemm_(CHARACTER TRANSA,CHARACTER TRANSB,INTEGER M,INTEGER N,INTEGER K,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZGEMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       CHARACTER TRANS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZGEMV  performs one of the matrix-vector operations<br>
*><br>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or<br>
*><br>
*>    y := alpha*A**H*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are vectors and A is an<br>
*> m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.<br>
*><br>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.<br>
*><br>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of DIMENSION at least<br>
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'<br>
*>           and at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.<br>
*>           Before entry with BETA non-zero, the incremented array Y<br>
*>           must contain the vector y. On exit, Y is overwritten by the<br>
*>           updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zgemv_(CHARACTER TRANS,INTEGER M,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZGERC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZGERC  performs the rank 1 operation<br>
*><br>
*>    A := alpha*x*y**H + A,<br>
*><br>
*> where alpha is a scalar, x is an m element vector, y is an n element<br>
*> vector and A is an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the m<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients. On exit, A is<br>
*>           overwritten by the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zgerc_(INTEGER M,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] A,INTEGER LDA);
/**
*> \brief \b ZGERU<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       INTEGER INCX,INCY,LDA,M,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZGERU  performs the rank 1 operation<br>
*><br>
*>    A := alpha*x*y**T + A,<br>
*><br>
*> where alpha is a scalar, x is an m element vector, y is an n element<br>
*> vector and A is an m by n matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of the matrix A.<br>
*>           M must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( m - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the m<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry, the leading m by n part of the array A must<br>
*>           contain the matrix of coefficients. On exit, A is<br>
*>           overwritten by the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zgeru_(INTEGER M,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] A,INTEGER LDA);
/**
*> \brief \b ZHBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER INCX,INCY,K,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n hermitian band matrix, with k super-diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the band matrix A is being supplied as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   The upper triangular part of A is<br>
*>                                  being supplied.<br>
*><br>
*>              UPLO = 'L' or 'l'   The lower triangular part of A is<br>
*>                                  being supplied.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry, K specifies the number of super-diagonals of the<br>
*>           matrix A. K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the hermitian matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer the upper<br>
*>           triangular part of a hermitian band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the hermitian matrix, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer the lower<br>
*>           triangular part of a hermitian band matrix from conventional<br>
*>           full matrix storage to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the<br>
*>           vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of DIMENSION at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the<br>
*>           vector y. On exit, Y is overwritten by the updated vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhbmv_(CHARACTER UPLO,INTEGER N,INTEGER K,double[] ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZHEMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER LDA,LDB,LDC,M,N<br>
*       CHARACTER SIDE,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*A*B + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*B*A + beta*C,<br>
*><br>
*> where alpha and beta are scalars, A is an hermitian matrix and  B and<br>
*> C are m by n matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE  specifies whether  the  hermitian matrix  A<br>
*>           appears on the  left or right  in the  operation as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,<br>
*><br>
*>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of  the  hermitian  matrix   A  is  to  be<br>
*>           referenced as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of the<br>
*>                                  hermitian matrix is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of the<br>
*>                                  hermitian matrix is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies the number of rows of the matrix  C.<br>
*>           M  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix C.<br>
*>           N  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           m  when  SIDE = 'L' or 'l'  and is n  otherwise.<br>
*>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of<br>
*>           the array  A  must contain the  hermitian matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading m by m upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  hermitian matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  m by m  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  hermitian<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of<br>
*>           the array  A  must contain the  hermitian matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading n by n upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  hermitian matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  n by n  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  hermitian<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Note that the imaginary parts  of the diagonal elements need<br>
*>           not be set, they are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array of DIMENSION ( LDB, n ).<br>
*>           Before entry, the leading  m by n part of the array  B  must<br>
*>           contain the matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n updated<br>
*>           matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhemm_(CHARACTER SIDE,CHARACTER UPLO,INTEGER M,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZHEMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEMV  performs the matrix-vector  operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n hermitian matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           lower triangular part of A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           upper triangular part of A is not referenced.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhemv_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZHER<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHER(UPLO,N,ALPHA,X,INCX,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHER   performs the hermitian rank 1 operation<br>
*><br>
*>    A := alpha*x*x**H + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
*> n by n hermitian matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           upper triangular part of A is not referenced. On exit, the<br>
*>           lower triangular part of the array A is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set, they are assumed to be zero, and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zher_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] A,INTEGER LDA);
/**
*> \brief \b ZHER2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       INTEGER INCX,INCY,LDA,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHER2  performs the hermitian rank 2 operation<br>
*><br>
*>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an n<br>
*> by n hermitian matrix.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           lower triangular part of A is not referenced. On exit, the<br>
*>           upper triangular part of the array A is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular part of the hermitian matrix and the strictly<br>
*>           upper triangular part of A is not referenced. On exit, the<br>
*>           lower triangular part of the array A is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set, they are assumed to be zero, and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zher2_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] A,INTEGER LDA);
/**
*> \brief \b ZHER2K<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       DOUBLE PRECISION BETA<br>
*       INTEGER K,LDA,LDB,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHER2K  performs one of the hermitian rank 2k operations<br>
*><br>
*>    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars with  beta  real,  C is an  n by n<br>
*> hermitian matrix and  A and B  are  n by k matrices in the first case<br>
*> and  k by n  matrices in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'    C := alpha*A*B**H          +<br>
*>                                         conjg( alpha )*B*A**H +<br>
*>                                         beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'    C := alpha*A**H*B          +<br>
*>                                         conjg( alpha )*B**H*A +<br>
*>                                         beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns  of the  matrices  A and B,  and on  entry  with<br>
*>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the<br>
*>           matrices  A and B.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16 .<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  k by n  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must<br>
*>           be at least  max( 1, k ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION .<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set,  they are assumed to be zero,  and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*><br>
*>  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.<br>
*>     Ed Anderson, Cray Research Inc.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zher2k_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZHERK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA,BETA<br>
*       INTEGER K,LDA,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHERK  performs one of the hermitian rank k operations<br>
*><br>
*>    C := alpha*A*A**H + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**H*A + beta*C,<br>
*><br>
*> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian<br>
*> matrix and  A  is an  n by k  matrix in the  first case and a  k by n<br>
*> matrix in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C.<br>
*><br>
*>              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the<br>
*>           matrix A.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION .<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION.<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  hermitian matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set,  they are assumed to be zero,  and on exit they<br>
*>           are set to zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*><br>
*>  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.<br>
*>     Ed Anderson, Cray Research Inc.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zherk_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,DOUBLE ALPHA,double[] A,INTEGER LDA,DOUBLE BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZHPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPMV  performs the matrix-vector operation<br>
*><br>
*>    y := alpha*A*x + beta*y,<br>
*><br>
*> where alpha and beta are scalars, x and y are n element vectors and<br>
*> A is an n by n hermitian matrix, supplied in packed form.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )<br>
*>           and a( 3, 1 ) respectively, and so on.<br>
*>           Note that the imaginary parts of the diagonal elements need<br>
*>           not be set and are assumed to be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta. When BETA is<br>
*>           supplied as zero then Y need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
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
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhpmv_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] AP,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
/**
*> \brief \b ZHPR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION ALPHA<br>
*       INTEGER INCX,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPR    performs the hermitian rank 1 operation<br>
*><br>
*>    A := alpha*x*x**H + A,<br>
*><br>
*> where alpha is a real scalar, x is an n element vector and A is an<br>
*> n by n hermitian matrix, supplied in packed form.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION.<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the hermitian matrix<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhpr_(CHARACTER UPLO,INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,double[] AP);
/**
*> \brief \b ZHPR2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       INTEGER INCX,INCY,N<br>
*       CHARACTER UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 AP(*),X(*),Y(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPR2  performs the hermitian rank 2 operation<br>
*><br>
*>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,<br>
*><br>
*> where alpha is a scalar, x and y are n element vectors and A is an<br>
*> n by n hermitian matrix, supplied in packed form.<br>
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
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] Y<br>
*> \verbatim<br>
*>          Y is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCY ) ).<br>
*>           Before entry, the incremented array Y must contain the n<br>
*>           element vector y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>           On entry, INCY specifies the increment for the elements of<br>
*>           Y. INCY must not be zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular part of the hermitian matrix<br>
*>           packed sequentially, column by column, so that AP( 1 )<br>
*>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )<br>
*>           and a( 2, 2 ) respectively, and so on. On exit, the array<br>
*>           AP is overwritten by the upper triangular part of the<br>
*>           updated matrix.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular part of the hermitian matrix<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhpr2_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] AP);
/**
*> \brief \b ZROTG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZROTG(CA,CB,C,S)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 CA,CB,S<br>
*       DOUBLE PRECISION C<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZROTG determines a double complex Givens rotation.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*  =====================================================================<br>
*/
	public void zrotg_(double[] CA,double[] CB,DOUBLE C,double[] S);
/**
*> \brief \b ZSCAL<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSCAL(N,ZA,ZX,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ZA<br>
*       INTEGER INCX,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZSCAL scales a vector by a constant.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 3/93 to return if incx .le. 0.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zscal_(INTEGER N,double[] ZA,double[] ZX,INTEGER INCX);
/**
*> \brief \b ZSWAP<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,INCY,N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 ZX(*),ZY(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZSWAP interchanges two vectors.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level1<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>     jack dongarra, 3/11/78.<br>
*>     modified 12/3/93, array(1) declarations changed to array(*)<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zswap_(INTEGER N,double[] ZX,INTEGER INCX,double[] ZY,INTEGER INCY);
/**
*> \brief \b ZSYMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER LDA,LDB,LDC,M,N<br>
*       CHARACTER SIDE,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYMM  performs one of the matrix-matrix operations<br>
*><br>
*>    C := alpha*A*B + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*B*A + beta*C,<br>
*><br>
*> where  alpha and beta are scalars, A is a symmetric matrix and  B and<br>
*> C are m by n matrices.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE  specifies whether  the  symmetric matrix  A<br>
*>           appears on the  left or right  in the  operation as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,<br>
*><br>
*>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of  the  symmetric  matrix   A  is  to  be<br>
*>           referenced as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of the<br>
*>                                  symmetric matrix is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry,  M  specifies the number of rows of the matrix  C.<br>
*>           M  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of the matrix C.<br>
*>           N  must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           m  when  SIDE = 'L' or 'l'  and is n  otherwise.<br>
*>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading m by m upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  m by m  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of<br>
*>           the array  A  must contain the  symmetric matrix,  such that<br>
*>           when  UPLO = 'U' or 'u', the leading n by n upper triangular<br>
*>           part of the array  A  must contain the upper triangular part<br>
*>           of the  symmetric matrix and the  strictly  lower triangular<br>
*>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',<br>
*>           the leading  n by n  lower triangular part  of the  array  A<br>
*>           must  contain  the  lower triangular part  of the  symmetric<br>
*>           matrix and the  strictly upper triangular part of  A  is not<br>
*>           referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then<br>
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at<br>
*>           least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array of DIMENSION ( LDB, n ).<br>
*>           Before entry, the leading  m by n part of the array  B  must<br>
*>           contain the matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is<br>
*>           supplied as zero then C need not be set on input.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry, the leading  m by n  part of the array  C must<br>
*>           contain the matrix  C,  except when  beta  is zero, in which<br>
*>           case C need not be set on entry.<br>
*>           On exit, the array  C  is overwritten by the  m by n updated<br>
*>           matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zsymm_(CHARACTER SIDE,CHARACTER UPLO,INTEGER M,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZSYR2K<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER K,LDA,LDB,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYR2K  performs one of the symmetric rank 2k operations<br>
*><br>
*>    C := alpha*A*B**T + alpha*B*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*B + alpha*B**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix<br>
*> and  A and B  are  n by k  matrices  in the  first  case  and  k by n<br>
*> matrices in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'    C := alpha*A*B**T + alpha*B*A**T +<br>
*>                                         beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'    C := alpha*A**T*B + alpha*B**T*A +<br>
*>                                         beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns  of the  matrices  A and B,  and on  entry  with<br>
*>           TRANS = 'T' or 't',  K  specifies  the number of rows of the<br>
*>           matrices  A and B.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  B  must contain the matrix  B,  otherwise<br>
*>           the leading  k by n  part of the array  B  must contain  the<br>
*>           matrix B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDB must be at least  max( 1, n ), otherwise  LDB must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zsyr2k_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZSYRK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA,BETA<br>
*       INTEGER K,LDA,LDC,N<br>
*       CHARACTER TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),C(LDC,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZSYRK  performs one of the symmetric rank k operations<br>
*><br>
*>    C := alpha*A*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*A + beta*C,<br>
*><br>
*> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix<br>
*> and  A  is an  n by k  matrix in the first case and a  k by n  matrix<br>
*> in the second case.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry,   UPLO  specifies  whether  the  upper  or  lower<br>
*>           triangular  part  of the  array  C  is to be  referenced  as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C<br>
*>                                  is to be referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry,  TRANS  specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'T' or 't',  K  specifies  the number of rows of the<br>
*>           matrix A.  K must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is<br>
*>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.<br>
*>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k<br>
*>           part of the array  A  must contain the matrix  A,  otherwise<br>
*>           the leading  k by n  part of the array  A  must contain  the<br>
*>           matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16<br>
*>           On entry, BETA specifies the scalar beta.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array of DIMENSION ( LDC, n ).<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n<br>
*>           upper triangular part of the array C must contain the upper<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           lower triangular part of C is not referenced.  On exit, the<br>
*>           upper triangular part of the array  C is overwritten by the<br>
*>           upper triangular part of the updated matrix.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n<br>
*>           lower triangular part of the array C must contain the lower<br>
*>           triangular part  of the  symmetric matrix  and the strictly<br>
*>           upper triangular part of C is not referenced.  On exit, the<br>
*>           lower triangular part of the array  C is overwritten by the<br>
*>           lower triangular part of the updated matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDC<br>
*> \verbatim<br>
*>          LDC is INTEGER<br>
*>           On entry, LDC specifies the first dimension of C as declared<br>
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zsyrk_(CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,double[] ALPHA,double[] A,INTEGER LDA,double[] BETA,double[] C,INTEGER LDC);
/**
*> \brief \b ZTBMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTBMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular band matrix, with ( k + 1 ) diagonals.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**H*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is (input/output) COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztbmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b ZTBSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,K,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTBSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )<br>
*> diagonals.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**H*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with UPLO = 'U' or 'u', K specifies the number of<br>
*>           super-diagonals of the matrix A.<br>
*>           On entry with UPLO = 'L' or 'l', K specifies the number of<br>
*>           sub-diagonals of the matrix A.<br>
*>           K must satisfy  0 .le. K.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the upper triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row<br>
*>           ( k + 1 ) of the array, the first super-diagonal starting at<br>
*>           position 2 in row k, and so on. The top left k by k triangle<br>
*>           of the array A is not referenced.<br>
*>           The following program segment will transfer an upper<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = K + 1 - J<br>
*>                    DO 10, I = MAX( 1, J - K ), J<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )<br>
*>           by n part of the array A must contain the lower triangular<br>
*>           band part of the matrix of coefficients, supplied column by<br>
*>           column, with the leading diagonal of the matrix in row 1 of<br>
*>           the array, the first sub-diagonal starting at position 1 in<br>
*>           row 2, and so on. The bottom right k by k triangle of the<br>
*>           array A is not referenced.<br>
*>           The following program segment will transfer a lower<br>
*>           triangular band matrix from conventional full matrix storage<br>
*>           to band storage:<br>
*><br>
*>                 DO 20, J = 1, N<br>
*>                    M = 1 - J<br>
*>                    DO 10, I = J, MIN( N, J + K )<br>
*>                       A( M + I, J ) = matrix( I, J )<br>
*>              10    CONTINUE<br>
*>              20 CONTINUE<br>
*><br>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A<br>
*>           corresponding to the diagonal elements of the matrix are not<br>
*>           referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           ( k + 1 ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztbsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER K,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b ZTPMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix, supplied in packed form.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**H*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is (input/output) COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztpmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] AP,double[] X,INTEGER INCX);
/**
*> \brief \b ZTPSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 AP(*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTPSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix, supplied in packed form.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**H*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array of DIMENSION at least<br>
*>           ( ( n*( n + 1 ) )/2 ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the array AP must<br>
*>           contain the upper triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )<br>
*>           respectively, and so on.<br>
*>           Before entry with UPLO = 'L' or 'l', the array AP must<br>
*>           contain the lower triangular matrix packed sequentially,<br>
*>           column by column, so that AP( 1 ) contains a( 1, 1 ),<br>
*>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )<br>
*>           respectively, and so on.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztpsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] AP,double[] X,INTEGER INCX);
/**
*> \brief \b ZTRMM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRMM  performs one of the matrix-matrix operations<br>
*><br>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )<br>
*><br>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from<br>
*>           the left or right as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, k ), where k is m<br>
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is (input/output) COMPLEX*16 array of DIMENSION ( LDB, n ).<br>
*>           Before entry,  the leading  m by n part of the array  B must<br>
*>           contain the matrix  B,  and  on exit  is overwritten  by the<br>
*>           transformed matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>           On entry, LDB specifies the first dimension of B as declared<br>
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least<br>
*>           max( 1, m ).<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrmm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB);
/**
*> \brief \b ZTRMV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRMV  performs one of the matrix-vector operations<br>
*><br>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,<br>
*><br>
*> where x is an n element vector and  A is an n by n unit, or non-unit,<br>
*> upper or lower triangular matrix.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   x := A*x.<br>
*><br>
*>              TRANS = 'T' or 't'   x := A**T*x.<br>
*><br>
*>              TRANS = 'C' or 'c'   x := A**H*x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in] X<br>
*> \verbatim<br>
*>          X is (input/output) COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element vector x. On exit, X is overwritten with the<br>
*>           tranformed vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrmv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,double[] X,INTEGER INCX);
/**
*> \brief \b ZTRSM<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       COMPLEX*16 ALPHA<br>
*       INTEGER LDA,LDB,M,N<br>
*       CHARACTER DIAG,SIDE,TRANSA,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),B(LDB,*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRSM  solves one of the matrix equations<br>
*><br>
*>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,<br>
*><br>
*> where alpha is a scalar, X and B are m by n matrices, A is a unit, or<br>
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of<br>
*><br>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.<br>
*><br>
*> The matrix X is overwritten on B.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>           On entry, SIDE specifies whether op( A ) appears on the left<br>
*>           or right of X as follows:<br>
*><br>
*>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.<br>
*><br>
*>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix A is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANSA<br>
*> \verbatim<br>
*>          TRANSA is CHARACTER*1<br>
*>           On entry, TRANSA specifies the form of op( A ) to be used in<br>
*>           the matrix multiplication as follows:<br>
*><br>
*>              TRANSA = 'N' or 'n'   op( A ) = A.<br>
*><br>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.<br>
*><br>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit triangular<br>
*>           as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>           On entry, M specifies the number of rows of B. M must be at<br>
*>           least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the number of columns of B.  N must be<br>
*>           at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16<br>
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is<br>
*>           zero then  A is not referenced and  B need not be set before<br>
*>           entry.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, k ),<br>
*>           where k is m when SIDE = 'L' or 'l'  <br>
*>             and k is n when SIDE = 'R' or 'r'.<br>
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k<br>
*>           upper triangular part of the array  A must contain the upper<br>
*>           triangular matrix  and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k<br>
*>           lower triangular part of the array  A must contain the lower<br>
*>           triangular matrix  and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of<br>
*>           A  are not referenced either,  but are assumed to be  unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then<br>
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'<br>
*>           then LDA must be at least max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is COMPLEX*16 array of DIMENSION ( LDB, n ).<br>
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
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level3<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 3 Blas routine.<br>
*><br>
*>  -- Written on 8-February-1989.<br>
*>     Jack Dongarra, Argonne National Laboratory.<br>
*>     Iain Duff, AERE Harwell.<br>
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.<br>
*>     Sven Hammarling, Numerical Algorithms Group Ltd.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrsm_(CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANSA,CHARACTER DIAG,INTEGER M,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] B,INTEGER LDB);
/**
*> \brief \b ZTRSV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER INCX,LDA,N<br>
*       CHARACTER DIAG,TRANS,UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16 A(LDA,*),X(*)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZTRSV  solves one of the systems of equations<br>
*><br>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,<br>
*><br>
*> where b and x are n element vectors and A is an n by n unit, or<br>
*> non-unit, upper or lower triangular matrix.<br>
*><br>
*> No test for singularity or near-singularity is included in this<br>
*> routine. Such tests must be performed before calling this routine.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On entry, UPLO specifies whether the matrix is an upper or<br>
*>           lower triangular matrix as follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.<br>
*><br>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the equations to be solved as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   A*x = b.<br>
*><br>
*>              TRANS = 'T' or 't'   A**T*x = b.<br>
*><br>
*>              TRANS = 'C' or 'c'   A**H*x = b.<br>
*> \endverbatim<br>
*><br>
*> \param[in] DIAG<br>
*> \verbatim<br>
*>          DIAG is CHARACTER*1<br>
*>           On entry, DIAG specifies whether or not A is unit<br>
*>           triangular as follows:<br>
*><br>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.<br>
*><br>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit<br>
*>                                  triangular.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix A.<br>
*>           N must be at least zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION ( LDA, n ).<br>
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n<br>
*>           upper triangular part of the array A must contain the upper<br>
*>           triangular matrix and the strictly lower triangular part of<br>
*>           A is not referenced.<br>
*>           Before entry with UPLO = 'L' or 'l', the leading n by n<br>
*>           lower triangular part of the array A must contain the lower<br>
*>           triangular matrix and the strictly upper triangular part of<br>
*>           A is not referenced.<br>
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of<br>
*>           A are not referenced either, but are assumed to be unity.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in the calling (sub) program. LDA must be at least<br>
*>           max( 1, n ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array of dimension at least<br>
*>           ( 1 + ( n - 1 )*abs( INCX ) ).<br>
*>           Before entry, the incremented array X must contain the n<br>
*>           element right-hand side vector b. On exit, X is overwritten<br>
*>           with the solution vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>           On entry, INCX specifies the increment for the elements of<br>
*>           X. INCX must not be zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16_blas_level2<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Level 2 Blas routine.<br>
*><br>
*>  -- Written on 22-October-1986.<br>
*>     Jack Dongarra, Argonne National Lab.<br>
*>     Jeremy Du Croz, Nag Central Office.<br>
*>     Sven Hammarling, Nag Central Office.<br>
*>     Richard Hanson, Sandia National Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ztrsv_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,double[] X,INTEGER INCX);

}