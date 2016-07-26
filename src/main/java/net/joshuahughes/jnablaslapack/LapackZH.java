package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZH extends Library
{

	public static LapackZH instance = (LapackZH) Native.loadLibrary("liblapack",LapackZH.class);

/**
*> \brief <b> ZHBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,<br>
*                         RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBEV computes all the eigenvalues and, optionally, eigenvectors of<br>
*> a complex Hermitian band matrix A.<br>
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
*>          AB is COMPLEX*16 array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, AB is overwritten by values generated during the<br>
*>          reduction to tridiagonal form.  If UPLO = 'U', the first<br>
*>          superdiagonal and the diagonal of the tridiagonal matrix T<br>
*>          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',<br>
*>          the diagonal and first subdiagonal of T are returned in the<br>
*>          first two rows of AB.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD + 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal<br>
*>          eigenvectors of the matrix A, with the i-th column of Z<br>
*>          holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (max(1,3*N-2))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of an intermediate tridiagonal<br>
*>                form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbev_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZHBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,<br>
*                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LRWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBEVD computes all the eigenvalues and, optionally, eigenvectors of<br>
*> a complex Hermitian band matrix A.  If eigenvectors are desired, it<br>
*> uses a divide and conquer algorithm.<br>
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
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
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
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, AB is overwritten by values generated during the<br>
*>          reduction to tridiagonal form.  If UPLO = 'U', the first<br>
*>          superdiagonal and the diagonal of the tridiagonal matrix T<br>
*>          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',<br>
*>          the diagonal and first subdiagonal of T are returned in the<br>
*>          first two rows of AB.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD + 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal<br>
*>          eigenvectors of the matrix A, with the i-th column of Z<br>
*>          holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
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
*>          If N <= 1,               LWORK must be at least 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK must be at least N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N**2.<br>
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
*>          The dimension of array RWORK.<br>
*>          If N <= 1,               LRWORK must be at least 1.<br>
*>          If JOBZ = 'N' and N > 1, LRWORK must be at least N.<br>
*>          If JOBZ = 'V' and N > 1, LRWORK must be at least<br>
*>                        1 + 5*N + 2*N**2.<br>
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
*>          The dimension of array IWORK.<br>
*>          If JOBZ = 'N' or N <= 1, LIWORK must be at least 1.<br>
*>          If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N .<br>
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
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of an intermediate tridiagonal<br>
*>                form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbevd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> ZHBEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL,<br>
*                          VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK,<br>
*                          IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a complex Hermitian band matrix A.  Eigenvalues and eigenvectors<br>
*> can be selected by specifying either a range of values or a range of<br>
*> indices for the desired eigenvalues.<br>
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
*>          = 'A': all eigenvalues will be found;<br>
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]<br>
*>                 will be found;<br>
*>          = 'I': the IL-th through IU-th eigenvalues will be found.<br>
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
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*><br>
*>          On exit, AB is overwritten by values generated during the<br>
*>          reduction to tridiagonal form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD + 1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDQ, N)<br>
*>          If JOBZ = 'V', the N-by-N unitary matrix used in the<br>
*>                          reduction to tridiagonal form.<br>
*>          If JOBZ = 'N', the array Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  If JOBZ = 'V', then<br>
*>          LDQ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing AB to tridiagonal form.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*DLAMCH('S').<br>
*><br>
*>          See "Computing Small Singular Values of Bidiagonal Matrices<br>
*>          with Guaranteed High Relative Accuracy," by Demmel and<br>
*>          Kahan, LAPACK Working Note #3.<br>
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
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If an eigenvector fails to converge, then that column of Z<br>
*>          contains the latest approximation to the eigenvector, and the<br>
*>          index of the eigenvector is returned in IFAIL.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvectors that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, then i eigenvectors failed to converge.<br>
*>                Their indices are stored in array IFAIL.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbevx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] Q,INTEGER LDQ,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZHBGST<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBGST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X,<br>
*                          LDX, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, VECT<br>
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBGST reduces a complex Hermitian-definite banded generalized<br>
*> eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,<br>
*> such that C has the same bandwidth as A.<br>
*><br>
*> B must have been previously factorized as S**H*S by ZPBSTF, using a<br>
*> split Cholesky factorization. A is overwritten by C = X**H*A*X, where<br>
*> X = S**(-1)*Q and Q is a unitary matrix chosen to preserve the<br>
*> bandwidth of A.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] VECT<br>
*> \verbatim<br>
*>          VECT is CHARACTER*1<br>
*>          = 'N':  do not form the transformation matrix X;<br>
*>          = 'V':  form X.<br>
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
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KA<br>
*> \verbatim<br>
*>          KA is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of superdiagonals of the matrix B if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KA >= KB >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first ka+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).<br>
*><br>
*>          On exit, the transformed matrix X**H*A*X, stored in the same<br>
*>          format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KA+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BB<br>
*> \verbatim<br>
*>          BB is COMPLEX*16 array, dimension (LDBB,N)<br>
*>          The banded factor S from the split Cholesky factorization of<br>
*>          B, as returned by ZPBSTF, stored in the first kb+1 rows of<br>
*>          the array.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDBB<br>
*> \verbatim<br>
*>          LDBB is INTEGER<br>
*>          The leading dimension of the array BB.  LDBB >= KB+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] X<br>
*> \verbatim<br>
*>          X is COMPLEX*16 array, dimension (LDX,N)<br>
*>          If VECT = 'V', the n-by-n matrix X.<br>
*>          If VECT = 'N', the array X is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDX<br>
*> \verbatim<br>
*>          LDX is INTEGER<br>
*>          The leading dimension of the array X.<br>
*>          LDX >= max(1,N) if VECT = 'V'; LDX >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
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
*> \ingroup complex16OTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbgst_(CHARACTER VECT,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,double[] AB,INTEGER LDAB,double[] BB,INTEGER LDBB,double[] X,INTEGER LDX,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHBGV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z,<br>
*                         LDZ, WORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBGV computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite banded eigenproblem, of<br>
*> the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian<br>
*> and banded, and B is also positive definite.<br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KA<br>
*> \verbatim<br>
*>          KA is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'. KA >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of superdiagonals of the matrix B if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'. KB >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first ka+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).<br>
*><br>
*>          On exit, the contents of AB are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KA+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BB<br>
*> \verbatim<br>
*>          BB is COMPLEX*16 array, dimension (LDBB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix B, stored in the first kb+1 rows of the array.  The<br>
*>          j-th column of B is stored in the j-th column of the array BB<br>
*>          as follows:<br>
*>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;<br>
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).<br>
*><br>
*>          On exit, the factor S from the split Cholesky factorization<br>
*>          B = S**H*S, as returned by ZPBSTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDBB<br>
*> \verbatim<br>
*>          LDBB is INTEGER<br>
*>          The leading dimension of the array BB.  LDBB >= KB+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors, with the i-th column of Z holding the<br>
*>          eigenvector associated with W(i). The eigenvectors are<br>
*>          normalized so that Z**H*B*Z = I.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is:<br>
*>             <= N:  the algorithm failed to converge:<br>
*>                    i off-diagonal elements of an intermediate<br>
*>                    tridiagonal form did not converge to zero;<br>
*>             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF<br>
*>                    returned INFO = i: B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbgv_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,double[] AB,INTEGER LDAB,double[] BB,INTEGER LDBB,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHBGVD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBGVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W,<br>
*                          Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK,<br>
*                          LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LRWORK,<br>
*      $                   LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBGVD computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite banded eigenproblem, of<br>
*> the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian<br>
*> and banded, and B is also positive definite.  If eigenvectors are<br>
*> desired, it uses a divide and conquer algorithm.<br>
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
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KA<br>
*> \verbatim<br>
*>          KA is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'. KA >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of superdiagonals of the matrix B if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'. KB >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first ka+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).<br>
*><br>
*>          On exit, the contents of AB are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KA+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BB<br>
*> \verbatim<br>
*>          BB is COMPLEX*16 array, dimension (LDBB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix B, stored in the first kb+1 rows of the array.  The<br>
*>          j-th column of B is stored in the j-th column of the array BB<br>
*>          as follows:<br>
*>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;<br>
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).<br>
*><br>
*>          On exit, the factor S from the split Cholesky factorization<br>
*>          B = S**H*S, as returned by ZPBSTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDBB<br>
*> \verbatim<br>
*>          LDBB is INTEGER<br>
*>          The leading dimension of the array BB.  LDBB >= KB+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors, with the i-th column of Z holding the<br>
*>          eigenvector associated with W(i). The eigenvectors are<br>
*>          normalized so that Z**H*B*Z = I.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO=0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If N <= 1,               LWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK >= N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK >= 2*N**2.<br>
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
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))<br>
*>          On exit, if INFO=0, RWORK(1) returns the optimal LRWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of array RWORK.<br>
*>          If N <= 1,               LRWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LRWORK >= N.<br>
*>          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.<br>
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
*>          On exit, if INFO=0, IWORK(1) returns the optimal LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of array IWORK.<br>
*>          If JOBZ = 'N' or N <= 1, LIWORK >= 1.<br>
*>          If JOBZ = 'V' and N > 1, LIWORK >= 3 + 5*N.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is:<br>
*>             <= N:  the algorithm failed to converge:<br>
*>                    i off-diagonal elements of an intermediate<br>
*>                    tridiagonal form did not converge to zero;<br>
*>             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF<br>
*>                    returned INFO = i: B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbgvd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,double[] AB,INTEGER LDAB,double[] BB,INTEGER LDBB,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZHBGVX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBGVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbgvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbgvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbgvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB,<br>
*                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,<br>
*                          LDZ, WORK, RWORK, IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M,<br>
*      $                   N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBGVX computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite banded eigenproblem, of<br>
*> the form A*x=(lambda)*B*x. Here A and B are assumed to be Hermitian<br>
*> and banded, and B is also positive definite.  Eigenvalues and<br>
*> eigenvectors can be selected by specifying either all eigenvalues,<br>
*> a range of values or a range of indices for the desired eigenvalues.<br>
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
*>          = 'A': all eigenvalues will be found;<br>
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]<br>
*>                 will be found;<br>
*>          = 'I': the IL-th through IU-th eigenvalues will be found.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KA<br>
*> \verbatim<br>
*>          KA is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'. KA >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of superdiagonals of the matrix B if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'. KB >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first ka+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).<br>
*><br>
*>          On exit, the contents of AB are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KA+1.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BB<br>
*> \verbatim<br>
*>          BB is COMPLEX*16 array, dimension (LDBB, N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix B, stored in the first kb+1 rows of the array.  The<br>
*>          j-th column of B is stored in the j-th column of the array BB<br>
*>          as follows:<br>
*>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;<br>
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).<br>
*><br>
*>          On exit, the factor S from the split Cholesky factorization<br>
*>          B = S**H*S, as returned by ZPBSTF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDBB<br>
*> \verbatim<br>
*>          LDBB is INTEGER<br>
*>          The leading dimension of the array BB.  LDBB >= KB+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDQ, N)<br>
*>          If JOBZ = 'V', the n-by-n matrix used in the reduction of<br>
*>          A*x = (lambda)*B*x to standard form, i.e. C*x = (lambda)*x,<br>
*>          and consequently C to tridiagonal form.<br>
*>          If JOBZ = 'N', the array Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  If JOBZ = 'N',<br>
*>          LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).<br>
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
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing AP to tridiagonal form.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*DLAMCH('S').<br>
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
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors, with the i-th column of Z holding the<br>
*>          eigenvector associated with W(i). The eigenvectors are<br>
*>          normalized so that Z**H*B*Z = I.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvectors that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, and i is:<br>
*>             <= N:  then i eigenvectors failed to converge.  Their<br>
*>                    indices are stored in array IFAIL.<br>
*>             > N:   if INFO = N + i, for 1 <= i <= N, then ZPBSTF<br>
*>                    returned INFO = i: B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void zhbgvx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,double[] AB,INTEGER LDAB,double[] BB,INTEGER LDBB,double[] Q,INTEGER LDQ,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZHBTRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHBTRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhbtrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhbtrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhbtrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, VECT<br>
*       INTEGER            INFO, KD, LDAB, LDQ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   D( * ), E( * )<br>
*       COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHBTRD reduces a complex Hermitian band matrix A to real symmetric<br>
*> tridiagonal form T by a unitary similarity transformation:<br>
*> Q**H * A * Q = T.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] VECT<br>
*> \verbatim<br>
*>          VECT is CHARACTER*1<br>
*>          = 'N':  do not form Q;<br>
*>          = 'V':  form Q;<br>
*>          = 'U':  update a matrix X, by forming X*Q.<br>
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
*> \param[in] KD<br>
*> \verbatim<br>
*>          KD is INTEGER<br>
*>          The number of superdiagonals of the matrix A if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is COMPLEX*16 array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the Hermitian band<br>
*>          matrix A, stored in the first KD+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).<br>
*>          On exit, the diagonal elements of AB are overwritten by the<br>
*>          diagonal elements of the tridiagonal matrix T; if KD > 0, the<br>
*>          elements on the first superdiagonal (if UPLO = 'U') or the<br>
*>          first subdiagonal (if UPLO = 'L') are overwritten by the<br>
*>          off-diagonal elements of T; the rest of AB is overwritten by<br>
*>          values generated during the reduction.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDAB<br>
*> \verbatim<br>
*>          LDAB is INTEGER<br>
*>          The leading dimension of the array AB.  LDAB >= KD+1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDQ,N)<br>
*>          On entry, if VECT = 'U', then Q must contain an N-by-N<br>
*>          matrix X; if VECT = 'N' or 'V', then Q need not be set.<br>
*><br>
*>          On exit:<br>
*>          if VECT = 'V', Q contains the N-by-N unitary matrix Q;<br>
*>          if VECT = 'U', Q contains the product X*Q;<br>
*>          if VECT = 'N', the array Q is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.<br>
*>          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.<br>
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
*> \date November 2011<br>
*<br>
*> \ingroup complex16OTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Modified by Linda Kaufman, Bell Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhbtrd_(CHARACTER VECT,CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] D,double[] E,double[] Q,INTEGER LDQ,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHECON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHECON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhecon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhecon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhecon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHECON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,<br>
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
*> ZHECON estimates the reciprocal of the condition number of a complex<br>
*> Hermitian matrix A using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by ZHETRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by ZHETRF.<br>
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
*>          as determined by ZHETRF.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhecon_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHECON_ROOK estimates the reciprocal of the condition number fort HE matrices using factorization obtained with one of the bounded diagonal pivoting methods (max 2 interchanges)<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHECON_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhecon_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhecon_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhecon_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHECON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,<br>
*                               INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       DOUBLE PRECISION   ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHECON_ROOK estimates the reciprocal of the condition number of a complex<br>
*> Hermitian matrix A using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by CHETRF_ROOK.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by CHETRF_ROOK.<br>
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
*>          as determined by CHETRF_ROOK.<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*> \verbatim<br>
*><br>
*>  November 2013,  Igor Kozachenko,<br>
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
	public void zhecon_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHEEQUB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEEQUB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheequb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheequb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheequb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )<br>
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
*> ZHEEQUB computes row and column scalings intended to equilibrate a<br>
*> Hermitian matrix A and reduce its condition number<br>
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
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
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
*>          The N-by-N Hermitian matrix whose scaling<br>
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
*> \date April 2012<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zheequb_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,double[] WORK,INTEGER INFO);
/**
*> \brief <b> ZHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,<br>
*                         INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEEV computes all eigenvalues and, optionally, eigenvectors of a<br>
*> complex Hermitian matrix A.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the<br>
*>          orthonormal eigenvectors of the matrix A.<br>
*>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')<br>
*>          or the upper triangle (if UPLO='U') of A, including the<br>
*>          diagonal, is destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          The length of the array WORK.  LWORK >= max(1,2*N-1).<br>
*>          For optimal efficiency, LWORK >= (NB+1)*N,<br>
*>          where NB is the blocksize for ZHETRD returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of an intermediate tridiagonal<br>
*>                form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zheev_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] W,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZHEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,<br>
*                          LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a<br>
*> complex Hermitian matrix A.  If eigenvectors are desired, it uses a<br>
*> divide and conquer algorithm.<br>
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
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the<br>
*>          orthonormal eigenvectors of the matrix A.<br>
*>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')<br>
*>          or the upper triangle (if UPLO='U') of A, including the<br>
*>          diagonal, is destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          The length of the array WORK.<br>
*>          If N <= 1,                LWORK must be at least 1.<br>
*>          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.<br>
*>          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.<br>
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
*>          If N <= 1,                LRWORK must be at least 1.<br>
*>          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.<br>
*>          If JOBZ  = 'V' and N > 1, LRWORK must be at least<br>
*>                         1 + 5*N + 2*N**2.<br>
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
*>          If N <= 1,                LIWORK must be at least 1.<br>
*>          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed<br>
*>                to converge; i off-diagonal elements of an intermediate<br>
*>                tridiagonal form did not converge to zero;<br>
*>                if INFO = i and JOBZ = 'V', then the algorithm failed<br>
*>                to compute an eigenvalue while working on the submatrix<br>
*>                lying in rows and columns INFO/(N+1) through<br>
*>                mod(INFO,N+1).<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*>  Modified description of INFO. Sven, 16 Feb 05.<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void zheevd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] W,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> ZHEEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEEVR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,<br>
*                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,<br>
*                          RWORK, LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK,<br>
*      $                   M, N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEEVR computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can<br>
*> be selected by specifying either a range of values or a range of<br>
*> indices for the desired eigenvalues.<br>
*><br>
*> ZHEEVR first reduces the matrix A to tridiagonal form T with a call<br>
*> to ZHETRD.  Then, whenever possible, ZHEEVR calls ZSTEMR to compute<br>
*> eigenspectrum using Relatively Robust Representations.  ZSTEMR<br>
*> computes eigenvalues by the dqds algorithm, while orthogonal<br>
*> eigenvectors are computed from various "good" L D L^T representations<br>
*> (also known as Relatively Robust Representations). Gram-Schmidt<br>
*> orthogonalization is avoided as far as possible. More specifically,<br>
*> the various steps of the algorithm are as follows.<br>
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
*> The desired accuracy of the output can be specified by the input<br>
*> parameter ABSTOL.<br>
*><br>
*> For more details, see DSTEMR's documentation and:<br>
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
*><br>
*> Note 1 : ZHEEVR calls ZSTEMR when the full spectrum is requested<br>
*> on machines which conform to the ieee-754 floating point standard.<br>
*> ZHEEVR calls DSTEBZ and ZSTEIN on non-ieee machines and<br>
*> when partial spectrum requests are made.<br>
*><br>
*> Normal execution of ZSTEMR may create NaNs and infinities and<br>
*> hence may abort due to a floating point exception in environments<br>
*> which do not handle NaNs and infinities in the ieee standard default<br>
*> manner.<br>
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
*>          For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and<br>
*>          ZSTEIN are called<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*>          On exit, the lower triangle (if UPLO='L') or the upper<br>
*>          triangle (if UPLO='U') of A, including the diagonal, is<br>
*>          destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing A to tridiagonal form.<br>
*><br>
*>          See "Computing Small Singular Values of Bidiagonal Matrices<br>
*>          with Guaranteed High Relative Accuracy," by Demmel and<br>
*>          Kahan, LAPACK Working Note #3.<br>
*><br>
*>          If high relative accuracy is important, set ABSTOL to<br>
*>          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that<br>
*>          eigenvalues are computed to high relative accuracy when<br>
*>          possible in future releases.  The current code does not<br>
*>          make any guarantees about high relative accuracy, but<br>
*>          furutre releases will. See J. Barlow and J. Demmel,<br>
*>          "Computing Accurate Eigensystems of Scaled Diagonally<br>
*>          Dominant Matrices", LAPACK Working Note #7, for a discussion<br>
*>          of which matrices define their eigenvalues to high relative<br>
*>          accuracy.<br>
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
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISUPPZ<br>
*> \verbatim<br>
*>          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) )<br>
*>          The support of the eigenvectors in Z, i.e., the indices<br>
*>          indicating the nonzero elements in Z. The i-th eigenvector<br>
*>          is nonzero only in elements ISUPPZ( 2*i-1 ) through<br>
*>          ISUPPZ( 2*i ).<br>
*>          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1<br>
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
*>          The length of the array WORK.  LWORK >= max(1,2*N).<br>
*>          For optimal efficiency, LWORK >= (NB+1)*N,<br>
*>          where NB is the max of the blocksize for ZHETRD and for<br>
*>          ZUNMTR as returned by ILAENV.<br>
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
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))<br>
*>          On exit, if INFO = 0, RWORK(1) returns the optimal<br>
*>          (and minimal) LRWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The length of the array RWORK.  LRWORK >= max(1,24*N).<br>
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
*>          On exit, if INFO = 0, IWORK(1) returns the optimal<br>
*>          (and minimal) LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.  LIWORK >= max(1,10*N).<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  Internal error<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Inderjit Dhillon, IBM Almaden, USA \n<br>
*>     Osni Marques, LBNL/NERSC, USA \n<br>
*>     Ken Stanley, Computer Science Division, University of<br>
*>       California at Berkeley, USA \n<br>
*>     Jason Riedy, Computer Science Division, University of<br>
*>       California at Berkeley, USA \n<br>
*><br>
*  =====================================================================<br>
*/
	public void zheevr_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,int[] ISUPPZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> ZHEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,<br>
*                          ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK,<br>
*                          IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can<br>
*> be selected by specifying either a range of values or a range of<br>
*> indices for the desired eigenvalues.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*>          On exit, the lower triangle (if UPLO='L') or the upper<br>
*>          triangle (if UPLO='U') of A, including the diagonal, is<br>
*>          destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing A to tridiagonal form.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*DLAMCH('S').<br>
*><br>
*>          See "Computing Small Singular Values of Bidiagonal Matrices<br>
*>          with Guaranteed High Relative Accuracy," by Demmel and<br>
*>          Kahan, LAPACK Working Note #3.<br>
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
*>          On normal exit, the first M elements contain the selected<br>
*>          eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If an eigenvector fails to converge, then that column of Z<br>
*>          contains the latest approximation to the eigenvector, and the<br>
*>          index of the eigenvector is returned in IFAIL.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
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
*>          The length of the array WORK.  LWORK >= 1, when N <= 1;<br>
*>          otherwise 2*N.<br>
*>          For optimal efficiency, LWORK >= (NB+1)*N,<br>
*>          where NB is the max of the blocksize for ZHETRD and for<br>
*>          ZUNMTR as returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvectors that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, then i eigenvectors failed to converge.<br>
*>                Their indices are stored in array IFAIL.<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zheevx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZHEGS2 reduces a Hermitian definite generalized eigenproblem to standard form, using the factorization results obtained from cpotrf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEGS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegs2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegs2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegs2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, N<br>
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
*> ZHEGS2 reduces a complex Hermitian-definite generalized<br>
*> eigenproblem to standard form.<br>
*><br>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,<br>
*> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)<br>
*><br>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or<br>
*> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H *A*L.<br>
*><br>
*> B must have been previously factorized as U**H *U or L*L**H by ZPOTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);<br>
*>          = 2 or 3: compute U*A*U**H or L**H *A*L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          Hermitian matrix A is stored, and how B has been factorized.<br>
*>          = 'U':  Upper triangular<br>
*>          = 'L':  Lower triangular<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          n by n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n by n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the transformed matrix, stored in the<br>
*>          same format as A.<br>
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
*>          The triangular factor from the Cholesky factorization of B,<br>
*>          as returned by ZPOTRF.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhegs2_(INTEGER ITYPE,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZHEGST<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEGST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, N<br>
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
*> ZHEGST reduces a complex Hermitian-definite generalized<br>
*> eigenproblem to standard form.<br>
*><br>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,<br>
*> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)<br>
*><br>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or<br>
*> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.<br>
*><br>
*> B must have been previously factorized as U**H*U or L*L**H by ZPOTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);<br>
*>          = 2 or 3: compute U*A*U**H or L**H*A*L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored and B is factored as<br>
*>                  U**H*U;<br>
*>          = 'L':  Lower triangle of A is stored and B is factored as<br>
*>                  L*L**H.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the transformed matrix, stored in the<br>
*>          same format as A.<br>
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
*>          The triangular factor from the Cholesky factorization of B,<br>
*>          as returned by ZPOTRF.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhegst_(INTEGER ITYPE,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZHEGV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,<br>
*                         LWORK, RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEGV computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.<br>
*> Here A and B are assumed to be Hermitian and B is also<br>
*> positive definite.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          Specifies the problem type to be solved:<br>
*>          = 1:  A*x = (lambda)*B*x<br>
*>          = 2:  A*B*x = (lambda)*x<br>
*>          = 3:  B*A*x = (lambda)*x<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*><br>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the<br>
*>          matrix Z of eigenvectors.  The eigenvectors are normalized<br>
*>          as follows:<br>
*>          if ITYPE = 1 or 2, Z**H*B*Z = I;<br>
*>          if ITYPE = 3, Z**H*inv(B)*Z = I.<br>
*>          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')<br>
*>          or the lower triangle (if UPLO='L') of A, including the<br>
*>          diagonal, is destroyed.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB, N)<br>
*>          On entry, the Hermitian positive definite matrix B.<br>
*>          If UPLO = 'U', the leading N-by-N upper triangular part of B<br>
*>          contains the upper triangular part of the matrix B.<br>
*>          If UPLO = 'L', the leading N-by-N lower triangular part of B<br>
*>          contains the lower triangular part of the matrix B.<br>
*><br>
*>          On exit, if INFO <= N, the part of B containing the matrix is<br>
*>          overwritten by the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**H*U or B = L*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          The length of the array WORK.  LWORK >= max(1,2*N-1).<br>
*>          For optimal efficiency, LWORK >= (NB+1)*N,<br>
*>          where NB is the blocksize for ZHETRD returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  ZPOTRF or ZHEEV returned an error code:<br>
*>             <= N:  if INFO = i, ZHEEV failed to converge;<br>
*>                    i off-diagonal elements of an intermediate<br>
*>                    tridiagonal form did not converge to zero;<br>
*>             > N:   if INFO = N + i, for 1 <= i <= N, then the leading<br>
*>                    minor of order i of B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhegv_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] W,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHEGVD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEGVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,<br>
*                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEGVD computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and<br>
*> B are assumed to be Hermitian and B is also positive definite.<br>
*> If eigenvectors are desired, it uses a divide and conquer algorithm.<br>
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
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          Specifies the problem type to be solved:<br>
*>          = 1:  A*x = (lambda)*B*x<br>
*>          = 2:  A*B*x = (lambda)*x<br>
*>          = 3:  B*A*x = (lambda)*x<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*><br>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the<br>
*>          matrix Z of eigenvectors.  The eigenvectors are normalized<br>
*>          as follows:<br>
*>          if ITYPE = 1 or 2, Z**H*B*Z = I;<br>
*>          if ITYPE = 3, Z**H*inv(B)*Z = I.<br>
*>          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')<br>
*>          or the lower triangle (if UPLO='L') of A, including the<br>
*>          diagonal, is destroyed.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB, N)<br>
*>          On entry, the Hermitian matrix B.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of B contains the<br>
*>          upper triangular part of the matrix B.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of B contains<br>
*>          the lower triangular part of the matrix B.<br>
*><br>
*>          On exit, if INFO <= N, the part of B containing the matrix is<br>
*>          overwritten by the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**H*U or B = L*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          The length of the array WORK.<br>
*>          If N <= 1,                LWORK >= 1.<br>
*>          If JOBZ  = 'N' and N > 1, LWORK >= N + 1.<br>
*>          If JOBZ  = 'V' and N > 1, LWORK >= 2*N + N**2.<br>
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
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))<br>
*>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of the array RWORK.<br>
*>          If N <= 1,                LRWORK >= 1.<br>
*>          If JOBZ  = 'N' and N > 1, LRWORK >= N.<br>
*>          If JOBZ  = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.<br>
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
*>          If N <= 1,                LIWORK >= 1.<br>
*>          If JOBZ  = 'N' and N > 1, LIWORK >= 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.<br>
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
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  ZPOTRF or ZHEEVD returned an error code:<br>
*>             <= N:  if INFO = i and JOBZ = 'N', then the algorithm<br>
*>                    failed to converge; i off-diagonal elements of an<br>
*>                    intermediate tridiagonal form did not converge to<br>
*>                    zero;<br>
*>                    if INFO = i and JOBZ = 'V', then the algorithm<br>
*>                    failed to compute an eigenvalue while working on<br>
*>                    the submatrix lying in rows and columns INFO/(N+1)<br>
*>                    through mod(INFO,N+1);<br>
*>             > N:   if INFO = N + i, for 1 <= i <= N, then the leading<br>
*>                    minor of order i of B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Modified so that no backsubstitution is performed if ZHEEVD fails to<br>
*>  converge (NEIG in old code could be greater than N causing out of<br>
*>  bounds reference to A - reported by Ralf Meyer).  Also corrected the<br>
*>  description of INFO and the test on ITYPE. Sven, 16 Feb 05.<br>
*> \endverbatim<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*><br>
*  =====================================================================<br>
*/
	public void zhegvd_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] W,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZHEGVX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHEGVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,<br>
*                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,<br>
*                          LWORK, RWORK, IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHEGVX computes selected eigenvalues, and optionally, eigenvectors<br>
*> of a complex generalized Hermitian-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and<br>
*> B are assumed to be Hermitian and B is also positive definite.<br>
*> Eigenvalues and eigenvectors can be selected by specifying either a<br>
*> range of values or a range of indices for the desired eigenvalues.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          Specifies the problem type to be solved:<br>
*>          = 1:  A*x = (lambda)*B*x<br>
*>          = 2:  A*B*x = (lambda)*x<br>
*>          = 3:  B*A*x = (lambda)*x<br>
*> \endverbatim<br>
*><br>
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
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
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
*>          A is COMPLEX*16 array, dimension (LDA, N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*><br>
*>          On exit,  the lower triangle (if UPLO='L') or the upper<br>
*>          triangle (if UPLO='U') of A, including the diagonal, is<br>
*>          destroyed.<br>
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
*>          B is COMPLEX*16 array, dimension (LDB, N)<br>
*>          On entry, the Hermitian matrix B.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of B contains the<br>
*>          upper triangular part of the matrix B.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of B contains<br>
*>          the lower triangular part of the matrix B.<br>
*><br>
*>          On exit, if INFO <= N, the part of B containing the matrix is<br>
*>          overwritten by the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**H*U or B = L*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDB<br>
*> \verbatim<br>
*>          LDB is INTEGER<br>
*>          The leading dimension of the array B.  LDB >= max(1,N).<br>
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
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing C to tridiagonal form, where C is the symmetric<br>
*>          matrix of the standard symmetric problem to which the<br>
*>          generalized problem is transformed.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*DLAMCH('S').<br>
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
*>          The first M elements contain the selected<br>
*>          eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          The eigenvectors are normalized as follows:<br>
*>          if ITYPE = 1 or 2, Z**T*B*Z = I;<br>
*>          if ITYPE = 3, Z**T*inv(B)*Z = I.<br>
*><br>
*>          If an eigenvector fails to converge, then that column of Z<br>
*>          contains the latest approximation to the eigenvector, and the<br>
*>          index of the eigenvector is returned in IFAIL.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
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
*>          The length of the array WORK.  LWORK >= max(1,2*N).<br>
*>          For optimal efficiency, LWORK >= (NB+1)*N,<br>
*>          where NB is the blocksize for ZHETRD returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvectors that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  ZPOTRF or ZHEEVX returned an error code:<br>
*>             <= N:  if INFO = i, ZHEEVX failed to converge;<br>
*>                    i eigenvectors failed to converge.  Their indices<br>
*>                    are stored in array IFAIL.<br>
*>             > N:   if INFO = N + i, for 1 <= i <= N, then the leading<br>
*>                    minor of order i of B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16HEeigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void zhegvx_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZHERFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHERFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zherfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zherfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zherfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,<br>
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
*> ZHERFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is Hermitian indefinite, and<br>
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
*>          AF is COMPLEX*16 array, dimension (LDAF,N)<br>
*>          The factored form of the matrix A.  AF contains the block<br>
*>          diagonal matrix D and the multipliers used to obtain the<br>
*>          factor U or L from the factorization A = U*D*U**H or<br>
*>          A = L*D*L**H as computed by ZHETRF.<br>
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
*>          as determined by ZHETRF.<br>
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
*>          On entry, the solution matrix X, as computed by ZHETRS.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zherfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHERFSX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHERFSX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zherfsx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zherfsx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zherfsx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHERFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
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
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZHERFSX improves the computed solution to a system of linear<br>
*>    equations when the coefficient matrix is Hermitian indefinite, and<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zherfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZHESV computes the solution to system of linear equations A * X = B for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHESV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,<br>
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
*> ZHESV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS<br>
*> matrices.<br>
*><br>
*> The diagonal pivoting method is used to factor A as<br>
*>    A = U * D * U**H,  if UPLO = 'U', or<br>
*>    A = L * D * L**H,  if UPLO = 'L',<br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is Hermitian and block diagonal with<br>
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
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the block diagonal matrix D and the<br>
*>          multipliers used to obtain the factor U or L from the<br>
*>          factorization A = U*D*U**H or A = L*D*L**H as computed by<br>
*>          ZHETRF.<br>
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
*>          determined by ZHETRF.  If IPIV(k) > 0, then rows and columns<br>
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
*>          ZHETRF.<br>
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
*> \ingroup complex16HEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zhesv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> ZHESVX computes the solution to system of linear equations A * X = B for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHESVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,<br>
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
*> ZHESVX uses the diagonal pivoting factorization to compute the<br>
*> solution to a complex system of linear equations A * X = B,<br>
*> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS<br>
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
*>       A = U * D * U**H,  if UPLO = 'U', or<br>
*>       A = L * D * L**H,  if UPLO = 'L',<br>
*>    where U (or L) is a product of permutation and unit upper (lower)<br>
*>    triangular matrices, and D is Hermitian and block diagonal with<br>
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
*> \param[in,out] AF<br>
*> \verbatim<br>
*>          AF is COMPLEX*16 array, dimension (LDAF,N)<br>
*>          If FACT = 'F', then AF is an input argument and on entry<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**H or A = L*D*L**H as computed by ZHETRF.<br>
*><br>
*>          If FACT = 'N', then AF is an output argument and on exit<br>
*>          returns the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**H or A = L*D*L**H.<br>
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
*>          of D, as determined by ZHETRF.<br>
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
*>          of D, as determined by ZHETRF.<br>
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
*>          NB is the optimal blocksize for ZHETRF.<br>
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
*> \ingroup complex16HEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zhesvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZHESVXX computes the solution to system of linear equations A * X = B for HE matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHESVXX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesvxx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesvxx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesvxx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHESVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
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
*      $                   WORK( * ), X( LDX, * )<br>
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
*>    ZHESVXX uses the diagonal pivoting factorization to compute the<br>
*>    solution to a complex*16 system of linear equations A * X = B, where<br>
*>    A is an N-by-N symmetric matrix and X and B are N-by-NRHS<br>
*>    matrices.<br>
*><br>
*>    If requested, both normwise and maximum componentwise error bounds<br>
*>    are returned. ZHESVXX will return a solution with a tiny<br>
*>    guaranteed error (O(eps) where eps is the working machine<br>
*>    precision) unless the matrix is very ill-conditioned, in which<br>
*>    case a warning is returned. Relevant condition numbers also are<br>
*>    calculated and returned.<br>
*><br>
*>    ZHESVXX accepts user-provided factorizations and equilibration<br>
*>    factors; see the definitions of the FACT and EQUED options.<br>
*>    Solving with refinement and using a factorization from a previous<br>
*>    ZHESVXX call will also produce a solution with either O(eps)<br>
*>    errors or warnings, but we cannot make that claim for general<br>
*>    user-provided factorizations and equilibration factors if they<br>
*>    differ from what ZHESVXX would itself produce.<br>
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
*>     structure of D, as determined by ZHETRF.  If IPIV(k) > 0,<br>
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
*>     structure of D, as determined by ZHETRF.<br>
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
*>          WORK is COMPLEX*16 array, dimension (5*N)<br>
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
*> \ingroup complex16HEsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void zhesvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE ZLA_HERPVGRW);
/**
*> \brief \b ZHESV_ROOK computes the solution to a system of linear equations A * X = B for HE matrices using the bounded Bunch-Kaufman ("rook") diagonal pivoting method<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHESV_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhesv_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhesv_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhesv_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHESV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,<br>
*                              LWORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHESV_ROOK computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS<br>
*> matrices.<br>
*><br>
*> The bounded Bunch-Kaufman ("rook") diagonal pivoting method is used<br>
*> to factor A as<br>
*>    A = U * D * U**T,  if UPLO = 'U', or<br>
*>    A = L * D * L**T,  if UPLO = 'L',<br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is Hermitian and block diagonal with<br>
*> 1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*> ZHETRF_ROOK is called to compute the factorization of a complex<br>
*> Hermition matrix A using the bounded Bunch-Kaufman ("rook") diagonal<br>
*> pivoting method.<br>
*><br>
*> The factored form of A is then used to solve the system<br>
*> of equations A * X = B by calling ZHETRS_ROOK (uses BLAS 2).<br>
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
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*><br>
*>          On exit, if INFO = 0, the block diagonal matrix D and the<br>
*>          multipliers used to obtain the factor U or L from the<br>
*>          factorization A = U*D*U**H or A = L*D*L**H as computed by<br>
*>          ZHETRF_ROOK.<br>
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
*>             Only the last KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             Only the first KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>             were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
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
*>          ZHETRF_ROOK.<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complex16HEsolve<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013,  Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*><br>
*>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,<br>
*>                  School of Mathematics,<br>
*>                  University of Manchester<br>
*><br>
*> \endverbatim<br>
*<br>
*<br>
*  =====================================================================<br>
*/
	public void zhesv_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZHESWAPR applies an elementary permutation on the rows and columns of a Hermitian matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHESWAPR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheswapr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheswapr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheswapr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHESWAPR( UPLO, N, A, LDA, I1, I2)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER        UPLO<br>
*       INTEGER          I1, I2, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16          A( LDA, N )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHESWAPR applies an elementary permutation on the rows and the columns of<br>
*> a hermitian matrix.<br>
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
*>          used to obtain the factor U or L as computed by CSYTRF.<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup complex16HEauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void zheswapr_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER I1,INTEGER I2);
/**
*> \brief \b ZHETD2 reduces a Hermitian matrix to real symmetric tridiagonal form by an unitary similarity transformation (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETD2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetd2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetd2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetd2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   D( * ), E( * )<br>
*       COMPLEX*16         A( LDA, * ), TAU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETD2 reduces a complex Hermitian matrix A to real symmetric<br>
*> tridiagonal form T by a unitary similarity transformation:<br>
*> Q**H * A * Q = T.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal<br>
*>          of A are overwritten by the corresponding elements of the<br>
*>          tridiagonal matrix T, and the elements above the first<br>
*>          superdiagonal, with the array TAU, represent the unitary<br>
*>          matrix Q as a product of elementary reflectors; if UPLO<br>
*>          = 'L', the diagonal and first subdiagonal of A are over-<br>
*>          written by the corresponding elements of the tridiagonal<br>
*>          matrix T, and the elements below the first subdiagonal, with<br>
*>          the array TAU, represent the unitary matrix Q as a product<br>
*>          of elementary reflectors. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX*16 array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \date September 2012<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(n-1) . . . H(2) H(1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in<br>
*>  A(1:i-1,i+1), and tau in TAU(i).<br>
*><br>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(n-1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),<br>
*>  and tau in TAU(i).<br>
*><br>
*>  The contents of A on exit are illustrated by the following examples<br>
*>  with n = 5:<br>
*><br>
*>  if UPLO = 'U':                       if UPLO = 'L':<br>
*><br>
*>    (  d   e   v2  v3  v4 )              (  d                  )<br>
*>    (      d   e   v3  v4 )              (  e   d              )<br>
*>    (          d   e   v4 )              (  v1  e   d          )<br>
*>    (              d   e  )              (  v1  v2  e   d      )<br>
*>    (                  d  )              (  v1  v2  v3  e   d  )<br>
*><br>
*>  where d and e denote diagonal and off-diagonal elements of T, and vi<br>
*>  denotes an element of the vector defining H(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhetd2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAU,INTEGER INFO);
/**
*> \brief \b ZHETF2 computes the factorization of a complex Hermitian matrix, using the diagonal pivoting method (unblocked algorithm, calling Level 2 BLAS).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHETF2 + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETF2( UPLO, N, A, LDA, IPIV, INFO )<br>
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
*> ZHETF2 computes the factorization of a complex Hermitian matrix A<br>
*> using the Bunch-Kaufman diagonal pivoting method:<br>
*><br>
*>    A = U*D*U**H  or  A = L*D*L**H<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, U**H is the conjugate transpose of U, and D is<br>
*> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**H, where<br>
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
*>  If UPLO = 'L', then A = L*D*L**H, where<br>
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
*>  09-29-06 - patch from<br>
*>    Bobby Cheng, MathWorks<br>
*><br>
*>    Replace l.210 and l.393<br>
*>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN<br>
*>    by<br>
*>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. DISNAN(ABSAKK) ) THEN<br>
*><br>
*>  01-01-96 - Based on modifications by<br>
*>    J. Lewis, Boeing Computer Services Company<br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b ZHETF2_ROOK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHETF2_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetf2_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetf2_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetf2_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETF2_ROOK( UPLO, N, A, LDA, IPIV, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16            A( LDA, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETF2_ROOK computes the factorization of a complex Hermitian matrix A<br>
*> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:<br>
*><br>
*>    A = U*D*U**H  or  A = L*D*L**H<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, U**H is the conjugate transpose of U, and D is<br>
*> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.<br>
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
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**H, where<br>
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
*>  If UPLO = 'L', then A = L*D*L**H, where<br>
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
*>  November 2013,  Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*><br>
*>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,<br>
*>                  School of Mathematics,<br>
*>                  University of Manchester<br>
*><br>
*>  01-01-96 - Based on modifications by<br>
*>    J. Lewis, Boeing Computer Services Company<br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetf2_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b ZHETRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   D( * ), E( * )<br>
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETRD reduces a complex Hermitian matrix A to real symmetric<br>
*> tridiagonal form T by a unitary similarity transformation:<br>
*> Q**H * A * Q = T.<br>
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
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal<br>
*>          of A are overwritten by the corresponding elements of the<br>
*>          tridiagonal matrix T, and the elements above the first<br>
*>          superdiagonal, with the array TAU, represent the unitary<br>
*>          matrix Q as a product of elementary reflectors; if UPLO<br>
*>          = 'L', the diagonal and first subdiagonal of A are over-<br>
*>          written by the corresponding elements of the tridiagonal<br>
*>          matrix T, and the elements below the first subdiagonal, with<br>
*>          the array TAU, represent the unitary matrix Q as a product<br>
*>          of elementary reflectors. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A.  LDA >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX*16 array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>          The dimension of the array WORK.  LWORK >= 1.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(n-1) . . . H(2) H(1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in<br>
*>  A(1:i-1,i+1), and tau in TAU(i).<br>
*><br>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(n-1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),<br>
*>  and tau in TAU(i).<br>
*><br>
*>  The contents of A on exit are illustrated by the following examples<br>
*>  with n = 5:<br>
*><br>
*>  if UPLO = 'U':                       if UPLO = 'L':<br>
*><br>
*>    (  d   e   v2  v3  v4 )              (  d                  )<br>
*>    (      d   e   v3  v4 )              (  e   d              )<br>
*>    (          d   e   v4 )              (  v1  e   d          )<br>
*>    (              d   e  )              (  v1  v2  e   d      )<br>
*>    (                  d  )              (  v1  v2  v3  e   d  )<br>
*><br>
*>  where d and e denote diagonal and off-diagonal elements of T, and vi<br>
*>  denotes an element of the vector defining H(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhetrd_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZHETRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
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
*> ZHETRF computes the factorization of a complex Hermitian matrix A<br>
*> using the Bunch-Kaufman diagonal pivoting method.  The form of the<br>
*> factorization is<br>
*><br>
*>    A = U*D*U**H  or  A = L*D*L**H<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is Hermitian and block diagonal with<br>
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
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  If UPLO = 'U', then A = U*D*U**H, where<br>
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
*>  If UPLO = 'L', then A = L*D*L**H, where<br>
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
	public void zhetrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZHETRF_ROOK computes the factorization of a complex Hermitian indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (blocked algorithm, calling Level 3 BLAS).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHETRF_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrf_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrf_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrf_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETRF_ROOK computes the factorization of a complex Hermitian matrix A<br>
*> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.<br>
*> The form of the factorization is<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is Hermitian and block diagonal with<br>
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
*>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading<br>
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
*>             Only the last KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were<br>
*>             interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k-1 and -IPIV(k-1) were inerchaged,<br>
*>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.<br>
*><br>
*>          If UPLO = 'L':<br>
*>             Only the first KB elements of IPIV are set.<br>
*><br>
*>             If IPIV(k) > 0, then rows and columns k and IPIV(k)<br>
*>             were interchanged and D(k,k) is a 1-by-1 diagonal block.<br>
*><br>
*>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and<br>
*>             columns k and -IPIV(k) were interchanged and rows and<br>
*>             columns k+1 and -IPIV(k+1) were inerchaged,<br>
*>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup complex16HEcomputational<br>
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
*>  June 2016,  Igor Kozachenko,<br>
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
	public void zhetrf_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZHETRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )<br>
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
*> ZHETRI computes the inverse of a complex Hermitian indefinite matrix<br>
*> A using the factorization A = U*D*U**H or A = L*D*L**H computed by<br>
*> ZHETRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          used to obtain the factor U or L as computed by ZHETRF.<br>
*><br>
*>          On exit, if INFO = 0, the (Hermitian) inverse of the original<br>
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
*>          as determined by ZHETRF.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetri_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHETRI2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRI2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
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
*> ZHETRI2 computes the inverse of a COMPLEX*16 hermitian indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> ZHETRF. ZHETRI2 set the LEADING DIMENSION of the workspace<br>
*> before calling ZHETRI2X that actually computes the inverse.<br>
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
*>          used to obtain the factor U or L as computed by ZHETRF.<br>
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
*>          as determined by ZHETRF.<br>
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
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>           calculates:<br>
*>              - the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array,<br>
*>              - and no error message related to LWORK is issued by XERBLA.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetri2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b ZHETRI2X<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRI2X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri2x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri2x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri2x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16            A( LDA, * ), WORK( N+NB+1,* )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETRI2X computes the inverse of a COMPLEX*16 Hermitian indefinite matrix<br>
*> A using the factorization A = U*D*U**H or A = L*D*L**H computed by<br>
*> ZHETRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          used to obtain the factor U or L as computed by ZHETRF.<br>
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
*>          as determined by ZHETRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N+NB+1,NB+3)<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetri2x_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER NB,INTEGER INFO);
/**
*> \brief \b ZHETRI_ROOK computes the inverse of HE matrix using the factorization obtained with the bounded Bunch-Kaufman ("rook") diagonal pivoting method.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHETRI_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetri_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetri_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetri_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX*16         A( LDA, * ), WORK( * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETRI_ROOK computes the inverse of a complex Hermitian indefinite matrix<br>
*> A using the factorization A = U*D*U**H or A = L*D*L**H computed by<br>
*> ZHETRF_ROOK.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          used to obtain the factor U or L as computed by ZHETRF_ROOK.<br>
*><br>
*>          On exit, if INFO = 0, the (Hermitian) inverse of the original<br>
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
*>          as determined by ZHETRF_ROOK.<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013,  Igor Kozachenko,<br>
*>                  Computer Science Division,<br>
*>                  University of California, Berkeley<br>
*><br>
*>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,<br>
*>                  School of Mathematics,<br>
*>                  University of Manchester<br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetri_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHETRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
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
*> ZHETRS solves a system of linear equations A*X = B with a complex<br>
*> Hermitian matrix A using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by ZHETRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by ZHETRF.<br>
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
*>          as determined by ZHETRF.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZHETRS2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHETRS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, <br>
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
*> ZHETRS2 solves a system of linear equations A*X = B with a complex<br>
*> Hermitian matrix A using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by ZHETRF and converted by ZSYCONV.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by ZHETRF.<br>
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
*>          as determined by ZHETRF.<br>
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
*> \ingroup complex16HEcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void zhetrs2_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHETRS_ROOK computes the solution to a system of linear equations A * X = B for HE matrices using factorization obtained with one of the bounded diagonal pivoting methods (max 2 interchanges)<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download ZHETRS_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrs_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrs_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrs_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       COMPLEX            A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHETRS_ROOK solves a system of linear equations A*X = B with a complex<br>
*> Hermitian matrix A using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by ZHETRF_ROOK.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by ZHETRF_ROOK.<br>
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
*>          as determined by ZHETRF_ROOK.<br>
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
*> \author Univ. of Tennessee<br>
*> \author Univ. of California Berkeley<br>
*> \author Univ. of Colorado Denver<br>
*> \author NAG Ltd.<br>
*<br>
*> \date November 2013<br>
*<br>
*> \ingroup complex16HEcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>  November 2013,  Igor Kozachenko,<br>
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
	public void zhetrs_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZHFRK performs a Hermitian rank-k operation for matrix in RFP format.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHFRK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhfrk.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhfrk.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhfrk.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA,<br>
*                         C )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION   ALPHA, BETA<br>
*       INTEGER            K, LDA, N<br>
*       CHARACTER          TRANS, TRANSR, UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * ), C( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> Level 3 BLAS like routine for C in RFP Format.<br>
*><br>
*> ZHFRK performs one of the Hermitian rank--k operations<br>
*><br>
*>    C := alpha*A*A**H + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**H*A + beta*C,<br>
*><br>
*> where alpha and beta are real scalars, C is an n--by--n Hermitian<br>
*> matrix and A is an n--by--k matrix in the first case and a k--by--n<br>
*> matrix in the second case.<br>
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
*><br>
*>           Unchanged on exit.<br>
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
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry,  N specifies the order of the matrix C.  N must be<br>
*>           at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number<br>
*>           of  columns   of  the   matrix   A,   and  on   entry   with<br>
*>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the<br>
*>           matrix A.  K must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is DOUBLE PRECISION<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array of DIMENSION (LDA,ka)<br>
*>           where KA<br>
*>           is K  when TRANS = 'N' or 'n', and is N otherwise. Before<br>
*>           entry with TRANS = 'N' or 'n', the leading N--by--K part of<br>
*>           the array A must contain the matrix A, otherwise the leading<br>
*>           K--by--N part of the array A must contain the matrix A.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>           On entry, LDA specifies the first dimension of A as declared<br>
*>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'<br>
*>           then  LDA must be at least  max( 1, n ), otherwise  LDA must<br>
*>           be at least  max( 1, k ).<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION<br>
*>           On entry, BETA specifies the scalar beta.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>           On entry, the matrix A in RFP Format. RFP Format is<br>
*>           described by TRANSR, UPLO and N. Note that the imaginary<br>
*>           parts of the diagonal elements need not be set, they are<br>
*>           assumed to be zero, and on exit they are set to zero.<br>
*> \endverbatim<br>
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
	public void zhfrk_(CHARACTER TRANSR,CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,DOUBLE ALPHA,double[] A,INTEGER LDA,DOUBLE BETA,double[] C);
/**
*> \brief \b ZHGEQZ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHGEQZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhgeqz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhgeqz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhgeqz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,<br>
*                          ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,<br>
*                          RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, COMPZ, JOB<br>
*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         ALPHA( * ), BETA( * ), H( LDH, * ),<br>
*      $                   Q( LDQ, * ), T( LDT, * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),<br>
*> where H is an upper Hessenberg matrix and T is upper triangular,<br>
*> using the single-shift QZ method.<br>
*> Matrix pairs of this type are produced by the reduction to<br>
*> generalized upper Hessenberg form of a complex matrix pair (A,B):<br>
*> <br>
*>    A = Q1*H*Z1**H,  B = Q1*T*Z1**H,<br>
*> <br>
*> as computed by ZGGHRD.<br>
*> <br>
*> If JOB='S', then the Hessenberg-triangular pair (H,T) is<br>
*> also reduced to generalized Schur form,<br>
*> <br>
*>    H = Q*S*Z**H,  T = Q*P*Z**H,<br>
*> <br>
*> where Q and Z are unitary matrices and S and P are upper triangular.<br>
*> <br>
*> Optionally, the unitary matrix Q from the generalized Schur<br>
*> factorization may be postmultiplied into an input matrix Q1, and the<br>
*> unitary matrix Z may be postmultiplied into an input matrix Z1.<br>
*> If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced<br>
*> the matrix pair (A,B) to generalized Hessenberg form, then the output<br>
*> matrices Q1*Q and Z1*Z are the unitary factors from the generalized<br>
*> Schur factorization of (A,B):<br>
*> <br>
*>    A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.<br>
*> <br>
*> To avoid overflow, eigenvalues of the matrix pair (H,T)<br>
*> (equivalently, of (A,B)) are computed as a pair of complex values<br>
*> (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an<br>
*> eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)<br>
*>    A*x = lambda*B*x<br>
*> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the<br>
*> alternate form of the GNEP<br>
*>    mu*A*y = B*y.<br>
*> The values of alpha and beta for the i-th eigenvalue can be read<br>
*> directly from the generalized Schur form:  alpha = S(i,i),<br>
*> beta = P(i,i).<br>
*><br>
*> Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix<br>
*>      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),<br>
*>      pp. 241--256.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          = 'E': Compute eigenvalues only;<br>
*>          = 'S': Computer eigenvalues and the Schur form.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          = 'N': Left Schur vectors (Q) are not computed;<br>
*>          = 'I': Q is initialized to the unit matrix and the matrix Q<br>
*>                 of left Schur vectors of (H,T) is returned;<br>
*>          = 'V': Q must contain a unitary matrix Q1 on entry and<br>
*>                 the product Q1*Q is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N': Right Schur vectors (Z) are not computed;<br>
*>          = 'I': Q is initialized to the unit matrix and the matrix Z<br>
*>                 of right Schur vectors of (H,T) is returned;<br>
*>          = 'V': Z must contain a unitary matrix Z1 on entry and<br>
*>                 the product Z1*Z is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices H, T, Q, and Z.  N >= 0.<br>
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
*>          ILO and IHI mark the rows and columns of H which are in<br>
*>          Hessenberg form.  It is assumed that A is already upper<br>
*>          triangular in rows and columns 1:ILO-1 and IHI+1:N.<br>
*>          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX*16 array, dimension (LDH, N)<br>
*>          On entry, the N-by-N upper Hessenberg matrix H.<br>
*>          On exit, if JOB = 'S', H contains the upper triangular<br>
*>          matrix S from the generalized Schur factorization.<br>
*>          If JOB = 'E', the diagonal of H matches that of S, but<br>
*>          the rest of H is unspecified.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H.  LDH >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] T<br>
*> \verbatim<br>
*>          T is COMPLEX*16 array, dimension (LDT, N)<br>
*>          On entry, the N-by-N upper triangular matrix T.<br>
*>          On exit, if JOB = 'S', T contains the upper triangular<br>
*>          matrix P from the generalized Schur factorization.<br>
*>          If JOB = 'E', the diagonal of T matches that of P, but<br>
*>          the rest of T is unspecified.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is COMPLEX*16 array, dimension (N)<br>
*>          The complex scalars alpha that define the eigenvalues of<br>
*>          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur<br>
*>          factorization.<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is COMPLEX*16 array, dimension (N)<br>
*>          The real non-negative scalars beta that define the<br>
*>          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized<br>
*>          Schur factorization.<br>
*><br>
*>          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)<br>
*>          represent the j-th eigenvalue of the matrix pair (A,B), in<br>
*>          one of the forms lambda = alpha/beta or mu = beta/alpha.<br>
*>          Since either lambda or mu may overflow, they should not,<br>
*>          in general, be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is COMPLEX*16 array, dimension (LDQ, N)<br>
*>          On entry, if COMPQ = 'V', the unitary matrix Q1 used in the<br>
*>          reduction of (A,B) to generalized Hessenberg form.<br>
*>          On exit, if COMPQ = 'I', the unitary matrix of left Schur<br>
*>          vectors of (H,T), and if COMPQ = 'V', the unitary matrix of<br>
*>          left Schur vectors of (A,B).<br>
*>          Not referenced if COMPQ = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= 1.<br>
*>          If COMPQ='V' or 'I', then LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the<br>
*>          reduction of (A,B) to generalized Hessenberg form.<br>
*>          On exit, if COMPZ = 'I', the unitary matrix of right Schur<br>
*>          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of<br>
*>          right Schur vectors of (A,B).<br>
*>          Not referenced if COMPZ = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1.<br>
*>          If COMPZ='V' or 'I', then LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.  LWORK >= max(1,N).<br>
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
*>          = 1,...,N: the QZ iteration did not converge.  (H,T) is not<br>
*>                     in Schur form, but ALPHA(i) and BETA(i),<br>
*>                     i=INFO+1,...,N should be correct.<br>
*>          = N+1,...,2*N: the shift calculation failed.  (H,T) is not<br>
*>                     in Schur form, but ALPHA(i) and BETA(i),<br>
*>                     i=INFO-N+1,...,N should be correct.<br>
*> \endverbatim<br>
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
*> \ingroup complex16GEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  We assume that complex ABS works as long as its value is less than<br>
*>  overflow.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhgeqz_(CHARACTER JOB,CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] T,INTEGER LDT,double[] ALPHA,double[] BETA,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHPCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )<br>
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
*> ZHPCON estimates the reciprocal of the condition number of a complex<br>
*> Hermitian packed matrix A using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by ZHPTRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by ZHPTRF, stored as a<br>
*>          packed triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZHPTRF.<br>
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
	public void zhpcon_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
/**
*> \brief <b> ZHPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,<br>
*                         INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPEV computes all the eigenvalues and, optionally, eigenvectors of a<br>
*> complex Hermitian matrix in packed storage.<br>
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
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, AP is overwritten by values generated during the<br>
*>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal<br>
*>          and first superdiagonal of the tridiagonal matrix T overwrite<br>
*>          the corresponding elements of A, and if UPLO = 'L', the<br>
*>          diagonal and first subdiagonal of T overwrite the<br>
*>          corresponding elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal<br>
*>          eigenvectors of the matrix A, with the i-th column of Z<br>
*>          holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (max(1, 2*N-1))<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of an intermediate tridiagonal<br>
*>                form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhpev_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] AP,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZHPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,<br>
*                          RWORK, LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDZ, LIWORK, LRWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPEVD computes all the eigenvalues and, optionally, eigenvectors of<br>
*> a complex Hermitian matrix A in packed storage.  If eigenvectors are<br>
*> desired, it uses a divide and conquer algorithm.<br>
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
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
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
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, AP is overwritten by values generated during the<br>
*>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal<br>
*>          and first superdiagonal of the tridiagonal matrix T overwrite<br>
*>          the corresponding elements of A, and if UPLO = 'L', the<br>
*>          diagonal and first subdiagonal of T overwrite the<br>
*>          corresponding elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal<br>
*>          eigenvectors of the matrix A, with the i-th column of Z<br>
*>          holding the eigenvector associated with W(i).<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the required LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of array WORK.<br>
*>          If N <= 1,               LWORK must be at least 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK must be at least N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the required sizes of the WORK, RWORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array,<br>
*>                                         dimension (LRWORK)<br>
*>          On exit, if INFO = 0, RWORK(1) returns the required LRWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of array RWORK.<br>
*>          If N <= 1,               LRWORK must be at least 1.<br>
*>          If JOBZ = 'N' and N > 1, LRWORK must be at least N.<br>
*>          If JOBZ = 'V' and N > 1, LRWORK must be at least<br>
*>                    1 + 5*N + 2*N**2.<br>
*><br>
*>          If LRWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the required sizes of the WORK, RWORK<br>
*>          and IWORK arrays, returns these values as the first entries<br>
*>          of the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of array IWORK.<br>
*>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the required sizes of the WORK, RWORK<br>
*>          and IWORK arrays, returns these values as the first entries<br>
*>          of the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of an intermediate tridiagonal<br>
*>                form did not converge to zero.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhpevd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] AP,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> ZHPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,<br>
*                          ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK,<br>
*                          IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, LDZ, M, N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a complex Hermitian matrix A in packed storage.<br>
*> Eigenvalues/vectors can be selected by specifying either a range of<br>
*> values or a range of indices for the desired eigenvalues.<br>
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
*>          = 'A': all eigenvalues will be found;<br>
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]<br>
*>                 will be found;<br>
*>          = 'I': the IL-th through IU-th eigenvalues will be found.<br>
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
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, AP is overwritten by values generated during the<br>
*>          reduction to tridiagonal form.  If UPLO = 'U', the diagonal<br>
*>          and first superdiagonal of the tridiagonal matrix T overwrite<br>
*>          the corresponding elements of A, and if UPLO = 'L', the<br>
*>          diagonal and first subdiagonal of T overwrite the<br>
*>          corresponding elements of A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is DOUBLE PRECISION<br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IL<br>
*> \verbatim<br>
*>          IL is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          smallest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing AP to tridiagonal form.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*DLAMCH('S').<br>
*><br>
*>          See "Computing Small Singular Values of Bidiagonal Matrices<br>
*>          with Guaranteed High Relative Accuracy," by Demmel and<br>
*>          Kahan, LAPACK Working Note #3.<br>
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
*>          If INFO = 0, the selected eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If an eigenvector fails to converge, then that column of Z<br>
*>          contains the latest approximation to the eigenvector, and<br>
*>          the index of the eigenvector is returned in IFAIL.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvectors that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, then i eigenvectors failed to converge.<br>
*>                Their indices are stored in array IFAIL.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhpevx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,double[] AP,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZHPGST<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPGST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITYPE, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         AP( * ), BP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPGST reduces a complex Hermitian-definite generalized<br>
*> eigenproblem to standard form, using packed storage.<br>
*><br>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,<br>
*> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)<br>
*><br>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or<br>
*> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.<br>
*><br>
*> B must have been previously factorized as U**H*U or L*L**H by ZPPTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);<br>
*>          = 2 or 3: compute U*A*U**H or L**H*A*L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored and B is factored as<br>
*>                  U**H*U;<br>
*>          = 'L':  Lower triangle of A is stored and B is factored as<br>
*>                  L*L**H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, if INFO = 0, the transformed matrix, stored in the<br>
*>          same format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] BP<br>
*> \verbatim<br>
*>          BP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          The triangular factor from the Cholesky factorization of B,<br>
*>          stored in the same format as A, as returned by ZPPTRF.<br>
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
	public void zhpgst_(INTEGER ITYPE,CHARACTER UPLO,INTEGER N,double[] AP,double[] BP,INTEGER INFO);
/**
*> \brief \b ZHPGV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,<br>
*                         RWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPGV computes all the eigenvalues and, optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.<br>
*> Here A and B are assumed to be Hermitian, stored in packed format,<br>
*> and B is also positive definite.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          Specifies the problem type to be solved:<br>
*>          = 1:  A*x = (lambda)*B*x<br>
*>          = 2:  A*B*x = (lambda)*x<br>
*>          = 3:  B*A*x = (lambda)*x<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the contents of AP are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BP<br>
*> \verbatim<br>
*>          BP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          B, packed columnwise in a linear array.  The j-th column of B<br>
*>          is stored in the array BP as follows:<br>
*>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**H*U or B = L*L**H, in the same storage<br>
*>          format as B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors.  The eigenvectors are normalized as follows:<br>
*>          if ITYPE = 1 or 2, Z**H*B*Z = I;<br>
*>          if ITYPE = 3, Z**H*inv(B)*Z = I.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (max(1, 2*N-1))<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  ZPPTRF or ZHPEV returned an error code:<br>
*>             <= N:  if INFO = i, ZHPEV failed to converge;<br>
*>                    i off-diagonal elements of an intermediate<br>
*>                    tridiagonal form did not convergeto zero;<br>
*>             > N:   if INFO = N + i, for 1 <= i <= n, then the leading<br>
*>                    minor of order i of B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void zhpgv_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] AP,double[] BP,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHPGVD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPGVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,<br>
*                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPGVD computes all the eigenvalues and, optionally, the eigenvectors<br>
*> of a complex generalized Hermitian-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and<br>
*> B are assumed to be Hermitian, stored in packed format, and B is also<br>
*> positive definite.<br>
*> If eigenvectors are desired, it uses a divide and conquer algorithm.<br>
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
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          Specifies the problem type to be solved:<br>
*>          = 1:  A*x = (lambda)*B*x<br>
*>          = 2:  A*B*x = (lambda)*x<br>
*>          = 3:  B*A*x = (lambda)*x<br>
*> \endverbatim<br>
*><br>
*> \param[in] JOBZ<br>
*> \verbatim<br>
*>          JOBZ is CHARACTER*1<br>
*>          = 'N':  Compute eigenvalues only;<br>
*>          = 'V':  Compute eigenvalues and eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the contents of AP are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BP<br>
*> \verbatim<br>
*>          BP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          B, packed columnwise in a linear array.  The j-th column of B<br>
*>          is stored in the array BP as follows:<br>
*>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**H*U or B = L*L**H, in the same storage<br>
*>          format as B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is DOUBLE PRECISION array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors.  The eigenvectors are normalized as follows:<br>
*>          if ITYPE = 1 or 2, Z**H*B*Z = I;<br>
*>          if ITYPE = 3, Z**H*inv(B)*Z = I.<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the required LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If N <= 1,               LWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK >= N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK >= 2*N.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the required sizes of the WORK, RWORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))<br>
*>          On exit, if INFO = 0, RWORK(1) returns the required LRWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LRWORK<br>
*> \verbatim<br>
*>          LRWORK is INTEGER<br>
*>          The dimension of array RWORK.<br>
*>          If N <= 1,               LRWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LRWORK >= N.<br>
*>          If JOBZ = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.<br>
*><br>
*>          If LRWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the required sizes of the WORK, RWORK<br>
*>          and IWORK arrays, returns these values as the first entries<br>
*>          of the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if INFO = 0, IWORK(1) returns the required LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of array IWORK.<br>
*>          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the required sizes of the WORK, RWORK<br>
*>          and IWORK arrays, returns these values as the first entries<br>
*>          of the WORK, RWORK and IWORK arrays, and no error message<br>
*>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  ZPPTRF or ZHPEVD returned an error code:<br>
*>             <= N:  if INFO = i, ZHPEVD failed to converge;<br>
*>                    i off-diagonal elements of an intermediate<br>
*>                    tridiagonal form did not convergeto zero;<br>
*>             > N:   if INFO = N + i, for 1 <= i <= n, then the leading<br>
*>                    minor of order i of B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void zhpgvd_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,double[] AP,double[] BP,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b ZHPGVX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPGVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,<br>
*                          IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK,<br>
*                          IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N<br>
*       DOUBLE PRECISION   ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       DOUBLE PRECISION   RWORK( * ), W( * )<br>
*       COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPGVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a complex generalized Hermitian-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and<br>
*> B are assumed to be Hermitian, stored in packed format, and B is also<br>
*> positive definite.  Eigenvalues and eigenvectors can be selected by<br>
*> specifying either a range of values or a range of indices for the<br>
*> desired eigenvalues.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          Specifies the problem type to be solved:<br>
*>          = 1:  A*x = (lambda)*B*x<br>
*>          = 2:  A*B*x = (lambda)*x<br>
*>          = 3:  B*A*x = (lambda)*x<br>
*> \endverbatim<br>
*><br>
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
*>          = 'A': all eigenvalues will be found;<br>
*>          = 'V': all eigenvalues in the half-open interval (VL,VU]<br>
*>                 will be found;<br>
*>          = 'I': the IL-th through IU-th eigenvalues will be found.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangles of A and B are stored;<br>
*>          = 'L':  Lower triangles of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices A and B.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the contents of AP are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] BP<br>
*> \verbatim<br>
*>          BP is COMPLEX*16 array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          B, packed columnwise in a linear array.  The j-th column of B<br>
*>          is stored in the array BP as follows:<br>
*>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**H*U or B = L*L**H, in the same storage<br>
*>          format as B.<br>
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
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IU<br>
*> \verbatim<br>
*>          IU is INTEGER<br>
*><br>
*>          If RANGE='I', the index of the<br>
*>          largest eigenvalue to be returned.<br>
*>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.<br>
*>          Not referenced if RANGE = 'A' or 'V'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ABSTOL<br>
*> \verbatim<br>
*>          ABSTOL is DOUBLE PRECISION<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less than<br>
*>          or equal to zero, then  EPS*|T|  will be used in its place,<br>
*>          where |T| is the 1-norm of the tridiagonal matrix obtained<br>
*>          by reducing AP to tridiagonal form.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*DLAMCH('S').<br>
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
*>          On normal exit, the first M elements contain the selected<br>
*>          eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ, N)<br>
*>          If JOBZ = 'N', then Z is not referenced.<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          The eigenvectors are normalized as follows:<br>
*>          if ITYPE = 1 or 2, Z**H*B*Z = I;<br>
*>          if ITYPE = 3, Z**H*inv(B)*Z = I.<br>
*><br>
*>          If an eigenvector fails to converge, then that column of Z<br>
*>          contains the latest approximation to the eigenvector, and the<br>
*>          index of the eigenvector is returned in IFAIL.<br>
*>          Note: the user must ensure that at least max(1,M) columns are<br>
*>          supplied in the array Z; if RANGE = 'V', the exact value of M<br>
*>          is not known in advance and an upper bound must be used.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1, and if<br>
*>          JOBZ = 'V', LDZ >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (2*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (7*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (N)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvectors that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  ZPPTRF or ZHPEVX returned an error code:<br>
*>             <= N:  if INFO = i, ZHPEVX failed to converge;<br>
*>                    i eigenvectors failed to converge.  Their indices<br>
*>                    are stored in array IFAIL.<br>
*>             > N:   if INFO = N + i, for 1 <= i <= n, then the leading<br>
*>                    minor of order i of B is not positive definite.<br>
*>                    The factorization of B could not be completed and<br>
*>                    no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup complex16OTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void zhpgvx_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,double[] AP,double[] BP,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,double[] WORK,double[] RWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b ZHPRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhprfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhprfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhprfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,<br>
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
*> ZHPRFS improves the computed solution to a system of linear<br>
*> equations when the coefficient matrix is Hermitian indefinite<br>
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
*>          The upper or lower triangle of the Hermitian matrix A, packed<br>
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
*>          factor U or L from the factorization A = U*D*U**H or<br>
*>          A = L*D*L**H as computed by ZHPTRF, stored as a packed<br>
*>          triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZHPTRF.<br>
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
*>          On entry, the solution matrix X, as computed by ZHPTRS.<br>
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
	public void zhprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief <b> ZHPSV computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )<br>
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
*> ZHPSV computes the solution to a complex system of linear equations<br>
*>    A * X = B,<br>
*> where A is an N-by-N Hermitian matrix stored in packed format and X<br>
*> and B are N-by-NRHS matrices.<br>
*><br>
*> The diagonal pivoting method is used to factor A as<br>
*>    A = U * D * U**H,  if UPLO = 'U', or<br>
*>    A = L * D * L**H,  if UPLO = 'L',<br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, D is Hermitian and block diagonal with 1-by-1<br>
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
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          See below for further details.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**H or A = L*D*L**H as computed by ZHPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D, as<br>
*>          determined by ZHPTRF.  If IPIV(k) > 0, then rows and columns<br>
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
	public void zhpsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> ZHPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X,<br>
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
*> ZHPSVX uses the diagonal pivoting factorization A = U*D*U**H or<br>
*> A = L*D*L**H to compute the solution to a complex system of linear<br>
*> equations A * X = B, where A is an N-by-N Hermitian matrix stored<br>
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
*>       A = U * D * U**H,  if UPLO = 'U', or<br>
*>       A = L * D * L**H,  if UPLO = 'L',<br>
*>    where U (or L) is a product of permutation and unit upper (lower)<br>
*>    triangular matrices and D is Hermitian and block diagonal with<br>
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
*>          = 'F':  On entry, AFP and IPIV contain the factored form of<br>
*>                  A.  AFP and IPIV will not be modified.<br>
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
*>          The upper or lower triangle of the Hermitian matrix A, packed<br>
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
*>          A = U*D*U**H or A = L*D*L**H as computed by ZHPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*><br>
*>          If FACT = 'N', then AFP is an output argument and on exit<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**H or A = L*D*L**H as computed by ZHPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          If FACT = 'F', then IPIV is an input argument and on entry<br>
*>          contains details of the interchanges and the block structure<br>
*>          of D, as determined by ZHPTRF.<br>
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
*>          of D, as determined by ZHPTRF.<br>
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
	public void zhpsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
/**
*> \brief \b ZHPTRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPTRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   D( * ), E( * )<br>
*       COMPLEX*16         AP( * ), TAU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHPTRD reduces a complex Hermitian matrix A stored in packed form to<br>
*> real symmetric tridiagonal form T by a unitary similarity<br>
*> transformation: Q**H * A * Q = T.<br>
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
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal<br>
*>          of A are overwritten by the corresponding elements of the<br>
*>          tridiagonal matrix T, and the elements above the first<br>
*>          superdiagonal, with the array TAU, represent the unitary<br>
*>          matrix Q as a product of elementary reflectors; if UPLO<br>
*>          = 'L', the diagonal and first subdiagonal of A are over-<br>
*>          written by the corresponding elements of the tridiagonal<br>
*>          matrix T, and the elements below the first subdiagonal, with<br>
*>          the array TAU, represent the unitary matrix Q as a product<br>
*>          of elementary reflectors. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is DOUBLE PRECISION array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is DOUBLE PRECISION array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is COMPLEX*16 array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*>  If UPLO = 'U', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(n-1) . . . H(2) H(1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,<br>
*>  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).<br>
*><br>
*>  If UPLO = 'L', the matrix Q is represented as a product of elementary<br>
*>  reflectors<br>
*><br>
*>     Q = H(1) H(2) . . . H(n-1).<br>
*><br>
*>  Each H(i) has the form<br>
*><br>
*>     H(i) = I - tau * v * v**H<br>
*><br>
*>  where tau is a complex scalar, and v is a complex vector with<br>
*>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,<br>
*>  overwriting A(i+2:n,i), and tau is stored in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhptrd_(CHARACTER UPLO,INTEGER N,double[] AP,double[] D,double[] E,double[] TAU,INTEGER INFO);
/**
*> \brief \b ZHPTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPTRF( UPLO, N, AP, IPIV, INFO )<br>
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
*> ZHPTRF computes the factorization of a complex Hermitian packed<br>
*> matrix A using the Bunch-Kaufman diagonal pivoting method:<br>
*><br>
*>    A = U*D*U**H  or  A = L*D*L**H<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is Hermitian and block diagonal with<br>
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
*>          On entry, the upper or lower triangle of the Hermitian matrix<br>
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
*>  If UPLO = 'U', then A = U*D*U**H, where<br>
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
*>  If UPLO = 'L', then A = L*D*L**H, where<br>
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
*>  J. Lewis, Boeing Computer Services Company<br>
*<br>
*  =====================================================================<br>
*/
	public void zhptrf_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,INTEGER INFO);
/**
*> \brief \b ZHPTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPTRI( UPLO, N, AP, IPIV, WORK, INFO )<br>
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
*> ZHPTRI computes the inverse of a complex Hermitian indefinite matrix<br>
*> A in packed storage using the factorization A = U*D*U**H or<br>
*> A = L*D*L**H computed by ZHPTRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          used to obtain the factor U or L as computed by ZHPTRF,<br>
*>          stored as a packed triangular matrix.<br>
*><br>
*>          On exit, if INFO = 0, the (Hermitian) inverse of the original<br>
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
*>          as determined by ZHPTRF.<br>
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
	public void zhptri_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,double[] WORK,INTEGER INFO);
/**
*> \brief \b ZHPTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHPTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhptrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhptrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhptrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )<br>
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
*> ZHPTRS solves a system of linear equations A*X = B with a complex<br>
*> Hermitian matrix A stored in packed format using the factorization<br>
*> A = U*D*U**H or A = L*D*L**H computed by ZHPTRF.<br>
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;<br>
*>          = 'L':  Lower triangular, form is A = L*D*L**H.<br>
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
*>          obtain the factor U or L as computed by ZHPTRF, stored as a<br>
*>          packed triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by ZHPTRF.<br>
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
	public void zhptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b ZHSEIN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHSEIN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhsein.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhsein.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhsein.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL,<br>
*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL,<br>
*                          IFAILR, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EIGSRC, INITV, SIDE<br>
*       INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IFAILL( * ), IFAILR( * )<br>
*       DOUBLE PRECISION   RWORK( * )<br>
*       COMPLEX*16         H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ZHSEIN uses inverse iteration to find specified right and/or left<br>
*> eigenvectors of a complex upper Hessenberg matrix H.<br>
*><br>
*> The right eigenvector x and the left eigenvector y of the matrix H<br>
*> corresponding to an eigenvalue w are defined by:<br>
*><br>
*>              H * x = w * x,     y**h * H = w * y**h<br>
*><br>
*> where y**h denotes the conjugate transpose of the vector y.<br>
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
*> \param[in] EIGSRC<br>
*> \verbatim<br>
*>          EIGSRC is CHARACTER*1<br>
*>          Specifies the source of eigenvalues supplied in W:<br>
*>          = 'Q': the eigenvalues were found using ZHSEQR; thus, if<br>
*>                 H has zero subdiagonal elements, and so is<br>
*>                 block-triangular, then the j-th eigenvalue can be<br>
*>                 assumed to be an eigenvalue of the block containing<br>
*>                 the j-th row/column.  This property allows ZHSEIN to<br>
*>                 perform inverse iteration on just one diagonal block.<br>
*>          = 'N': no assumptions are made on the correspondence<br>
*>                 between eigenvalues and diagonal blocks.  In this<br>
*>                 case, ZHSEIN must always perform inverse iteration<br>
*>                 using the whole matrix H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INITV<br>
*> \verbatim<br>
*>          INITV is CHARACTER*1<br>
*>          = 'N': no initial vectors are supplied;<br>
*>          = 'U': user-supplied initial vectors are stored in the arrays<br>
*>                 VL and/or VR.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          Specifies the eigenvectors to be computed. To select the<br>
*>          eigenvector corresponding to the eigenvalue W(j),<br>
*>          SELECT(j) must be set to .TRUE..<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] H<br>
*> \verbatim<br>
*>          H is COMPLEX*16 array, dimension (LDH,N)<br>
*>          The upper Hessenberg matrix H.<br>
*>          If a NaN is detected in H, the routine will return with INFO=-6.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H.  LDH >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] W<br>
*> \verbatim<br>
*>          W is COMPLEX*16 array, dimension (N)<br>
*>          On entry, the eigenvalues of H.<br>
*>          On exit, the real parts of W may have been altered since<br>
*>          close eigenvalues are perturbed slightly in searching for<br>
*>          independent eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)<br>
*>          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must<br>
*>          contain starting vectors for the inverse iteration for the<br>
*>          left eigenvectors; the starting vector for each eigenvector<br>
*>          must be in the same column in which the eigenvector will be<br>
*>          stored.<br>
*>          On exit, if SIDE = 'L' or 'B', the left eigenvectors<br>
*>          specified by SELECT will be stored consecutively in the<br>
*>          columns of VL, in the same order as their eigenvalues.<br>
*>          If SIDE = 'R', VL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.<br>
*>          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)<br>
*>          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must<br>
*>          contain starting vectors for the inverse iteration for the<br>
*>          right eigenvectors; the starting vector for each eigenvector<br>
*>          must be in the same column in which the eigenvector will be<br>
*>          stored.<br>
*>          On exit, if SIDE = 'R' or 'B', the right eigenvectors<br>
*>          specified by SELECT will be stored consecutively in the<br>
*>          columns of VR, in the same order as their eigenvalues.<br>
*>          If SIDE = 'L', VR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.<br>
*>          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.<br>
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
*>          The number of columns in the arrays VL and/or VR required to<br>
*>          store the eigenvectors (= the number of .TRUE. elements in<br>
*>          SELECT).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (N*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] RWORK<br>
*> \verbatim<br>
*>          RWORK is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAILL<br>
*> \verbatim<br>
*>          IFAILL is INTEGER array, dimension (MM)<br>
*>          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left<br>
*>          eigenvector in the i-th column of VL (corresponding to the<br>
*>          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the<br>
*>          eigenvector converged satisfactorily.<br>
*>          If SIDE = 'R', IFAILL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAILR<br>
*> \verbatim<br>
*>          IFAILR is INTEGER array, dimension (MM)<br>
*>          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right<br>
*>          eigenvector in the i-th column of VR (corresponding to the<br>
*>          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the<br>
*>          eigenvector converged satisfactorily.<br>
*>          If SIDE = 'L', IFAILR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, i is the number of eigenvectors which<br>
*>                failed to converge; see IFAILL and IFAILR for further<br>
*>                details.<br>
*> \endverbatim<br>
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
*>  Each eigenvector is normalized so that the element of largest<br>
*>  magnitude has magnitude 1; here the magnitude of a complex number<br>
*>  (x,y) is taken to be |x|+|y|.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void zhsein_(CHARACTER SIDE,CHARACTER EIGSRC,CHARACTER INITV,boolean[] SELECT,INTEGER N,double[] H,INTEGER LDH,double[] W,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,double[] RWORK,int[] IFAILL,int[] IFAILR,INTEGER INFO);
/**
*> \brief \b ZHSEQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZHSEQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhseqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhseqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhseqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,<br>
*                          WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N<br>
*       CHARACTER          COMPZ, JOB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    ZHSEQR computes the eigenvalues of a Hessenberg matrix H<br>
*>    and, optionally, the matrices T and Z from the Schur decomposition<br>
*>    H = Z T Z**H, where T is an upper triangular matrix (the<br>
*>    Schur form), and Z is the unitary matrix of Schur vectors.<br>
*><br>
*>    Optionally Z may be postmultiplied into an input unitary<br>
*>    matrix Q so that this routine can give the Schur factorization<br>
*>    of a matrix A which has been reduced to the Hessenberg form H<br>
*>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>           = 'E':  compute eigenvalues only;<br>
*>           = 'S':  compute eigenvalues and the Schur form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>           = 'N':  no Schur vectors are computed;<br>
*>           = 'I':  Z is initialized to the unit matrix and the matrix Z<br>
*>                   of Schur vectors of H is returned;<br>
*>           = 'V':  Z must contain an unitary matrix Q on entry, and<br>
*>                   the product Q*Z is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The order of the matrix H.  N .GE. 0.<br>
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
*>           It is assumed that H is already upper triangular in rows<br>
*>           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally<br>
*>           set by a previous call to ZGEBAL, and then passed to ZGEHRD<br>
*>           when the matrix output by ZGEBAL is reduced to Hessenberg<br>
*>           form. Otherwise ILO and IHI should be set to 1 and N<br>
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.<br>
*>           If N = 0, then ILO = 1 and IHI = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is COMPLEX*16 array, dimension (LDH,N)<br>
*>           On entry, the upper Hessenberg matrix H.<br>
*>           On exit, if INFO = 0 and JOB = 'S', H contains the upper<br>
*>           triangular matrix T from the Schur decomposition (the<br>
*>           Schur form). If INFO = 0 and JOB = 'E', the contents of<br>
*>           H are unspecified on exit.  (The output value of H when<br>
*>           INFO.GT.0 is given under the description of INFO below.)<br>
*><br>
*>           Unlike earlier versions of ZHSEQR, this subroutine may<br>
*>           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1<br>
*>           or j = IHI+1, IHI+2, ... N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>           The leading dimension of the array H. LDH .GE. max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is COMPLEX*16 array, dimension (N)<br>
*>           The computed eigenvalues. If JOB = 'S', the eigenvalues are<br>
*>           stored in the same order as on the diagonal of the Schur<br>
*>           form returned in H, with W(i) = H(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is COMPLEX*16 array, dimension (LDZ,N)<br>
*>           If COMPZ = 'N', Z is not referenced.<br>
*>           If COMPZ = 'I', on entry Z need not be set and on exit,<br>
*>           if INFO = 0, Z contains the unitary matrix Z of the Schur<br>
*>           vectors of H.  If COMPZ = 'V', on entry Z must contain an<br>
*>           N-by-N matrix Q, which is assumed to be equal to the unit<br>
*>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,<br>
*>           if INFO = 0, Z contains Q*Z.<br>
*>           Normally Q is the unitary matrix generated by ZUNGHR<br>
*>           after the call to ZGEHRD which formed the Hessenberg matrix<br>
*>           H. (The output value of Z when INFO.GT.0 is given under<br>
*>           the description of INFO below.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>           The leading dimension of the array Z.  if COMPZ = 'I' or<br>
*>           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is COMPLEX*16 array, dimension (LWORK)<br>
*>           On exit, if INFO = 0, WORK(1) returns an estimate of<br>
*>           the optimal value for LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)<br>
*>           is sufficient and delivers very good and sometimes<br>
*>           optimal performance.  However, LWORK as large as 11*N<br>
*>           may be required for optimal performance.  A workspace<br>
*>           query is recommended to determine the optimal workspace<br>
*>           size.<br>
*><br>
*>           If LWORK = -1, then ZHSEQR does a workspace query.<br>
*>           In this case, ZHSEQR checks the input parameters and<br>
*>           estimates the optimal workspace size for the given<br>
*>           values of N, ILO and IHI.  The estimate is returned<br>
*>           in WORK(1).  No error message related to LWORK is<br>
*>           issued by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>             =  0:  successful exit<br>
*>           .LT. 0:  if INFO = -i, the i-th argument had an illegal<br>
*>                    value<br>
*>           .GT. 0:  if INFO = i, ZHSEQR failed to compute all of<br>
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR<br>
*>                and WI contain those eigenvalues which have been<br>
*>                successfully computed.  (Failures are rare.)<br>
*><br>
*>                If INFO .GT. 0 and JOB = 'E', then on exit, the<br>
*>                remaining unconverged eigenvalues are the eigen-<br>
*>                values of the upper Hessenberg matrix rows and<br>
*>                columns ILO through INFO of the final, output<br>
*>                value of H.<br>
*><br>
*>                If INFO .GT. 0 and JOB   = 'S', then on exit<br>
*><br>
*>           (*)  (initial value of H)*U  = U*(final value of H)<br>
*><br>
*>                where U is a unitary matrix.  The final<br>
*>                value of  H is upper Hessenberg and triangular in<br>
*>                rows and columns INFO+1 through IHI.<br>
*><br>
*>                If INFO .GT. 0 and COMPZ = 'V', then on exit<br>
*><br>
*>                  (final value of Z)  =  (initial value of Z)*U<br>
*><br>
*>                where U is the unitary matrix in (*) (regard-<br>
*>                less of the value of JOB.)<br>
*><br>
*>                If INFO .GT. 0 and COMPZ = 'I', then on exit<br>
*>                      (final value of Z)  = U<br>
*>                where U is the unitary matrix in (*) (regard-<br>
*>                less of the value of JOB.)<br>
*><br>
*>                If INFO .GT. 0 and COMPZ = 'N', then Z is not<br>
*>                accessed.<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>             Default values supplied by<br>
*>             ILAENV(ISPEC,'ZHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).<br>
*>             It is suggested that these defaults be adjusted in order<br>
*>             to attain best performance in each particular<br>
*>             computational environment.<br>
*><br>
*>            ISPEC=12: The ZLAHQR vs ZLAQR0 crossover point.<br>
*>                      Default: 75. (Must be at least 11.)<br>
*><br>
*>            ISPEC=13: Recommended deflation window size.<br>
*>                      This depends on ILO, IHI and NS.  NS is the<br>
*>                      number of simultaneous shifts returned<br>
*>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)<br>
*>                      The default for (IHI-ILO+1).LE.500 is NS.<br>
*>                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.<br>
*><br>
*>            ISPEC=14: Nibble crossover point. (See IPARMQ for<br>
*>                      details.)  Default: 14% of deflation window<br>
*>                      size.<br>
*><br>
*>            ISPEC=15: Number of simultaneous shifts in a multishift<br>
*>                      QR iteration.<br>
*><br>
*>                      If IHI-ILO+1 is ...<br>
*><br>
*>                      greater than      ...but less    ... the<br>
*>                      or equal to ...      than        default is<br>
*><br>
*>                           1               30          NS =   2(+)<br>
*>                          30               60          NS =   4(+)<br>
*>                          60              150          NS =  10(+)<br>
*>                         150              590          NS =  **<br>
*>                         590             3000          NS =  64<br>
*>                        3000             6000          NS = 128<br>
*>                        6000             infinity      NS = 256<br>
*><br>
*>                  (+)  By default some or all matrices of this order<br>
*>                       are passed to the implicit double shift routine<br>
*>                       ZLAHQR and this parameter is ignored.  See<br>
*>                       ISPEC=12 above and comments in IPARMQ for<br>
*>                       details.<br>
*><br>
*>                 (**)  The asterisks (**) indicate an ad-hoc<br>
*>                       function of N increasing from 10 to 64.<br>
*><br>
*>            ISPEC=16: Select structured matrix multiply.<br>
*>                      If the number of simultaneous shifts (specified<br>
*>                      by ISPEC=15) is less than 14, then the default<br>
*>                      for ISPEC=16 is 0.  Otherwise the default for<br>
*>                      ISPEC=16 is 2.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3<br>
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages<br>
*>       929--947, 2002.<br>
*> \n<br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal<br>
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.<br>
*<br>
*  =====================================================================<br>
*/
	public void zhseqr_(CHARACTER JOB,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] W,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);

}