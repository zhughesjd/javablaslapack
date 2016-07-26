package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSS extends Library
{

	public static LapackSS instance = (LapackSS) Native.loadLibrary("liblapack",LapackSS.class);

/**
*> \brief <b> SSBEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,<br>
*                         INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBEV computes all the eigenvalues and, optionally, eigenvectors of<br>
*> a real symmetric band matrix A.<br>
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
*>          AB is REAL array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
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
*>          WORK is REAL array, dimension (max(1,3*N-2))<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbev_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief <b> SSBEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,<br>
*                          LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KD, LDAB, LDZ, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AB( LDAB, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBEVD computes all the eigenvalues and, optionally, eigenvectors of<br>
*> a real symmetric band matrix A. If eigenvectors are desired, it uses<br>
*> a divide and conquer algorithm.<br>
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
*>          AB is REAL array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
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
*>          WORK is REAL array,<br>
*>                                         dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          IF N <= 1,                LWORK must be at least 1.<br>
*>          If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.<br>
*>          If JOBZ  = 'V' and N > 2, LWORK must be at least<br>
*>                         ( 1 + 5*N + 2*N**2 ).<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
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
*>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.<br>
*>          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbevd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> SSBEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL,<br>
*                          VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK,<br>
*                          IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric band matrix A.  Eigenvalues and eigenvectors can<br>
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
*>          AB is REAL array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*> \param[out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ, N)<br>
*>          If JOBZ = 'V', the N-by-N orthogonal matrix used in the<br>
*>                         reduction to tridiagonal form.<br>
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
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M))<br>
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
*>          WORK is REAL array, dimension (7*N)<br>
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
*>          = 0:  successful exit.<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbevx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] Q,INTEGER LDQ,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSBGST<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBGST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbgst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbgst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbgst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X,<br>
*                          LDX, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, VECT<br>
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), BB( LDBB, * ), WORK( * ),<br>
*      $                   X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBGST reduces a real symmetric-definite banded generalized<br>
*> eigenproblem  A*x = lambda*B*x  to standard form  C*y = lambda*y,<br>
*> such that C has the same bandwidth as A.<br>
*><br>
*> B must have been previously factorized as S**T*S by SPBSTF, using a<br>
*> split Cholesky factorization. A is overwritten by C = X**T*A*X, where<br>
*> X = S**(-1)*Q and Q is an orthogonal matrix chosen to preserve the<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix A, stored in the first ka+1 rows of the array.  The<br>
*>          j-th column of A is stored in the j-th column of the array AB<br>
*>          as follows:<br>
*>          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;<br>
*>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).<br>
*><br>
*>          On exit, the transformed matrix X**T*A*X, stored in the same<br>
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
*>          BB is REAL array, dimension (LDBB,N)<br>
*>          The banded factor S from the split Cholesky factorization of<br>
*>          B, as returned by SPBSTF, stored in the first KB+1 rows of<br>
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
*>          X is REAL array, dimension (LDX,N)<br>
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
*>          WORK is REAL array, dimension (2*N)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbgst_(CHARACTER VECT,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,float[] AB,INTEGER LDAB,float[] BB,INTEGER LDBB,float[] X,INTEGER LDX,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSBGV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbgv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbgv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbgv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z,<br>
*                         LDZ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), BB( LDBB, * ), W( * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBGV computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a real generalized symmetric-definite banded eigenproblem, of<br>
*> the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric<br>
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
*>          AB is REAL array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*>          BB is REAL array, dimension (LDBB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix B, stored in the first kb+1 rows of the array.  The<br>
*>          j-th column of B is stored in the j-th column of the array BB<br>
*>          as follows:<br>
*>          if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;<br>
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).<br>
*><br>
*>          On exit, the factor S from the split Cholesky factorization<br>
*>          B = S**T*S, as returned by SPBSTF.<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors, with the i-th column of Z holding the<br>
*>          eigenvector associated with W(i). The eigenvectors are<br>
*>          normalized so that Z**T*B*Z = I.<br>
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
*>          WORK is REAL array, dimension (3*N)<br>
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
*>             > N:   if INFO = N + i, for 1 <= i <= N, then SPBSTF<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbgv_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,float[] AB,INTEGER LDAB,float[] BB,INTEGER LDBB,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSBGVD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBGVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbgvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbgvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbgvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W,<br>
*                          Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AB( LDAB, * ), BB( LDBB, * ), W( * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBGVD computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a real generalized symmetric-definite banded eigenproblem, of the<br>
*> form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and<br>
*> banded, and B is also positive definite.  If eigenvectors are<br>
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
*>          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] KB<br>
*> \verbatim<br>
*>          KB is INTEGER<br>
*>          The number of superdiagonals of the matrix B if UPLO = 'U',<br>
*>          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is REAL array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*>          BB is REAL array, dimension (LDBB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix B, stored in the first kb+1 rows of the array.  The<br>
*>          j-th column of B is stored in the j-th column of the array BB<br>
*>          as follows:<br>
*>          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;<br>
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).<br>
*><br>
*>          On exit, the factor S from the split Cholesky factorization<br>
*>          B = S**T*S, as returned by SPBSTF.<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors, with the i-th column of Z holding the<br>
*>          eigenvector associated with W(i).  The eigenvectors are<br>
*>          normalized so Z**T*B*Z = I.<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If N <= 1,               LWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK >= 3*N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.<br>
*>          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
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
*>             > N:   if INFO = N + i, for 1 <= i <= N, then SPBSTF<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbgvd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,float[] AB,INTEGER LDAB,float[] BB,INTEGER LDBB,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b SSBGVX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBGVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbgvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbgvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbgvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB,<br>
*                          LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,<br>
*                          LDZ, WORK, IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, KA, KB, LDAB, LDBB, LDQ, LDZ, M,<br>
*      $                   N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               AB( LDAB, * ), BB( LDBB, * ), Q( LDQ, * ),<br>
*      $                   W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBGVX computes selected eigenvalues, and optionally, eigenvectors<br>
*> of a real generalized symmetric-definite banded eigenproblem, of<br>
*> the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric<br>
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
*>          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AB<br>
*> \verbatim<br>
*>          AB is REAL array, dimension (LDAB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*>          BB is REAL array, dimension (LDBB, N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
*>          matrix B, stored in the first kb+1 rows of the array.  The<br>
*>          j-th column of B is stored in the j-th column of the array BB<br>
*>          as follows:<br>
*>          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;<br>
*>          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).<br>
*><br>
*>          On exit, the factor S from the split Cholesky factorization<br>
*>          B = S**T*S, as returned by SPBSTF.<br>
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
*>          Q is REAL array, dimension (LDQ, N)<br>
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
*>          VL is REAL<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors, with the i-th column of Z holding the<br>
*>          eigenvector associated with W(i).  The eigenvectors are<br>
*>          normalized so Z**T*B*Z = I.<br>
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
*>          WORK is REAL array, dimension (7N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (5N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAIL<br>
*> \verbatim<br>
*>          IFAIL is INTEGER array, dimension (M)<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M elements of<br>
*>          IFAIL are zero.  If INFO > 0, then IFAIL contains the<br>
*>          indices of the eigenvalues that failed to converge.<br>
*>          If JOBZ = 'N', then IFAIL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0 : successful exit<br>
*>          < 0 : if INFO = -i, the i-th argument had an illegal value<br>
*>          <= N: if INFO = i, then i eigenvectors failed to converge.<br>
*>                  Their indices are stored in IFAIL.<br>
*>          > N : SPBSTF returned an error code; i.e.,<br>
*>                if INFO = N + i, for 1 <= i <= N, then the leading<br>
*>                minor of order i of B is not positive definite.<br>
*>                The factorization of B could not be completed and<br>
*>                no eigenvalues or eigenvectors were computed.<br>
*> \endverbatim<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void ssbgvx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,INTEGER KA,INTEGER KB,float[] AB,INTEGER LDAB,float[] BB,INTEGER LDBB,float[] Q,INTEGER LDQ,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSBTRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSBTRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssbtrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssbtrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssbtrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,<br>
*                          WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, VECT<br>
*       INTEGER            INFO, KD, LDAB, LDQ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AB( LDAB, * ), D( * ), E( * ), Q( LDQ, * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSBTRD reduces a real symmetric band matrix A to symmetric<br>
*> tridiagonal form T by an orthogonal similarity transformation:<br>
*> Q**T * A * Q = T.<br>
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
*>          AB is REAL array, dimension (LDAB,N)<br>
*>          On entry, the upper or lower triangle of the symmetric band<br>
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
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is REAL array, dimension (LDQ,N)<br>
*>          On entry, if VECT = 'U', then Q must contain an N-by-N<br>
*>          matrix X; if VECT = 'N' or 'V', then Q need not be set.<br>
*><br>
*>          On exit:<br>
*>          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;<br>
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
*>          WORK is REAL array, dimension (N)<br>
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
*>  Modified by Linda Kaufman, Bell Labs.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssbtrd_(CHARACTER VECT,CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] D,float[] E,float[] Q,INTEGER LDQ,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSFRK performs a symmetric rank-k operation for matrix in RFP format.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSFRK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssfrk.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssfrk.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssfrk.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA,<br>
*                         C )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       REAL               ALPHA, BETA<br>
*       INTEGER            K, LDA, N<br>
*       CHARACTER          TRANS, TRANSR, UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), C( * )<br>
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
*> SSFRK performs one of the symmetric rank--k operations<br>
*><br>
*>    C := alpha*A*A**T + beta*C,<br>
*><br>
*> or<br>
*><br>
*>    C := alpha*A**T*A + beta*C,<br>
*><br>
*> where alpha and beta are real scalars, C is an n--by--n symmetric<br>
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
*>          = 'T':  The Transpose Form of RFP A is stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>           On  entry, UPLO specifies whether the upper or lower<br>
*>           triangular part of the array C is to be referenced as<br>
*>           follows:<br>
*><br>
*>              UPLO = 'U' or 'u'   Only the upper triangular part of C<br>
*>                                  is to be referenced.<br>
*><br>
*>              UPLO = 'L' or 'l'   Only the lower triangular part of C<br>
*>                                  is to be referenced.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] TRANS<br>
*> \verbatim<br>
*>          TRANS is CHARACTER*1<br>
*>           On entry, TRANS specifies the operation to be performed as<br>
*>           follows:<br>
*><br>
*>              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.<br>
*><br>
*>              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.<br>
*><br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           On entry, N specifies the order of the matrix C. N must be<br>
*>           at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] K<br>
*> \verbatim<br>
*>          K is INTEGER<br>
*>           On entry with TRANS = 'N' or 'n', K specifies the number<br>
*>           of  columns of the matrix A, and on entry with TRANS = 'T'<br>
*>           or 't', K specifies the number of rows of the matrix A. K<br>
*>           must be at least zero.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ALPHA<br>
*> \verbatim<br>
*>          ALPHA is REAL<br>
*>           On entry, ALPHA specifies the scalar alpha.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array of DIMENSION (LDA,ka)<br>
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
*>          BETA is REAL<br>
*>           On entry, BETA specifies the scalar beta.<br>
*>           Unchanged on exit.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] C<br>
*> \verbatim<br>
*>          C is REAL array, dimension (NT)<br>
*>           NT = N*(N+1)/2. On entry, the symmetric matrix C in RFP<br>
*>           Format. RFP Format is described by TRANSR, UPLO and N.<br>
*> \endverbatim<br>
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
	public void ssfrk_(CHARACTER TRANSR,CHARACTER UPLO,CHARACTER TRANS,INTEGER N,INTEGER K,REAL ALPHA,float[] A,INTEGER LDA,REAL BETA,float[] C);
/**
*> \brief \b SSPCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspcon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspcon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspcon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPCON estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a real symmetric packed matrix A using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by SSPTRF.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSPTRF, stored as a<br>
*>          packed triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by SSPTRF.<br>
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
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
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
	public void sspcon_(CHARACTER UPLO,INTEGER N,float[] AP,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> SSPEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPEV computes all the eigenvalues and, optionally, eigenvectors of a<br>
*> real symmetric matrix A in packed storage.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
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
*>          WORK is REAL array, dimension (3*N)<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sspev_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] AP,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief <b> SSPEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,<br>
*                          IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDZ, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPEVD computes all the eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric matrix A in packed storage. If eigenvectors are<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the required LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If N <= 1,               LWORK must be at least 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK must be at least<br>
*>                                                 1 + 6*N + N**2.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the required sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
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
*>          The dimension of the array IWORK.<br>
*>          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the required sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sspevd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] AP,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> SSPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,<br>
*                          ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, LDZ, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric matrix A in packed storage.  Eigenvalues/vectors<br>
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
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the selected eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M))<br>
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
*>          WORK is REAL array, dimension (8*N)<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sspevx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,float[] AP,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSPGST<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPGST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPGST( ITYPE, UPLO, N, AP, BP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITYPE, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), BP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPGST reduces a real symmetric-definite generalized eigenproblem<br>
*> to standard form, using packed storage.<br>
*><br>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,<br>
*> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)<br>
*><br>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or<br>
*> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.<br>
*><br>
*> B must have been previously factorized as U**T*U or L*L**T by SPPTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);<br>
*>          = 2 or 3: compute U*A*U**T or L**T*A*L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored and B is factored as<br>
*>                  U**T*U;<br>
*>          = 'L':  Lower triangle of A is stored and B is factored as<br>
*>                  L*L**T.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          BP is REAL array, dimension (N*(N+1)/2)<br>
*>          The triangular factor from the Cholesky factorization of B,<br>
*>          stored in the same format as A, as returned by SPPTRF.<br>
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
	public void sspgst_(INTEGER ITYPE,CHARACTER UPLO,INTEGER N,float[] AP,float[] BP,INTEGER INFO);
/**
*> \brief \b SSPGV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,<br>
*                         INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), BP( * ), W( * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPGV computes all the eigenvalues and, optionally, the eigenvectors<br>
*> of a real generalized symmetric-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.<br>
*> Here A and B are assumed to be symmetric, stored in packed format,<br>
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
*>          AP is REAL array, dimension<br>
*>                            (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          BP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          B, packed columnwise in a linear array.  The j-th column of B<br>
*>          is stored in the array BP as follows:<br>
*>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**T*U or B = L*L**T, in the same storage<br>
*>          format as B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors.  The eigenvectors are normalized as follows:<br>
*>          if ITYPE = 1 or 2, Z**T*B*Z = I;<br>
*>          if ITYPE = 3, Z**T*inv(B)*Z = I.<br>
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
*>          WORK is REAL array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  SPPTRF or SSPEV returned an error code:<br>
*>             <= N:  if INFO = i, SSPEV failed to converge;<br>
*>                    i off-diagonal elements of an intermediate<br>
*>                    tridiagonal form did not converge to zero.<br>
*>             > N:   if INFO = n + i, for 1 <= i <= n, then the leading<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sspgv_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] AP,float[] BP,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSPGVD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPGVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,<br>
*                          LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDZ, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               AP( * ), BP( * ), W( * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPGVD computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a real generalized symmetric-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and<br>
*> B are assumed to be symmetric, stored in packed format, and B is also<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          BP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          B, packed columnwise in a linear array.  The j-th column of B<br>
*>          is stored in the array BP as follows:<br>
*>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**T*U or B = L*L**T, in the same storage<br>
*>          format as B.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of<br>
*>          eigenvectors.  The eigenvectors are normalized as follows:<br>
*>          if ITYPE = 1 or 2, Z**T*B*Z = I;<br>
*>          if ITYPE = 3, Z**T*inv(B)*Z = I.<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the required LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If N <= 1,               LWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK >= 2*N.<br>
*>          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the required sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
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
*>          The dimension of the array IWORK.<br>
*>          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.<br>
*>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the required sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  SPPTRF or SSPEVD returned an error code:<br>
*>             <= N:  if INFO = i, SSPEVD failed to converge;<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void sspgvd_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] AP,float[] BP,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b SSPGVX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPGVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspgvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspgvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspgvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,<br>
*                          IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK,<br>
*                          IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               AP( * ), BP( * ), W( * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPGVX computes selected eigenvalues, and optionally, eigenvectors<br>
*> of a real generalized symmetric-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A<br>
*> and B are assumed to be symmetric, stored in packed storage, and B<br>
*> is also positive definite.  Eigenvalues and eigenvectors can be<br>
*> selected by specifying either a range of values or a range of indices<br>
*> for the desired eigenvalues.<br>
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
*>          = 'U':  Upper triangle of A and B are stored;<br>
*>          = 'L':  Lower triangle of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix pencil (A,B).  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] AP<br>
*> \verbatim<br>
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
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
*>          BP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          B, packed columnwise in a linear array.  The j-th column of B<br>
*>          is stored in the array BP as follows:<br>
*>          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.<br>
*><br>
*>          On exit, the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**T*U or B = L*L**T, in the same storage<br>
*>          format as B.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          On normal exit, the first M elements contain the selected<br>
*>          eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M))<br>
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
*>          WORK is REAL array, dimension (8*N)<br>
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
*>          > 0:  SPPTRF or SSPEVX returned an error code:<br>
*>             <= N:  if INFO = i, SSPEVX failed to converge;<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void sspgvx_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,float[] AP,float[] BP,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSPRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssprfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssprfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssprfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,<br>
*                          FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               AFP( * ), AP( * ), B( LDB, * ), BERR( * ),<br>
*      $                   FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPRFS improves the computed solution to a system of linear<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The upper or lower triangle of the symmetric matrix A, packed<br>
*>          columnwise in a linear array.  The j-th column of A is stored<br>
*>          in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*> \endverbatim<br>
*><br>
*> \param[in] AFP<br>
*> \verbatim<br>
*>          AFP is REAL array, dimension (N*(N+1)/2)<br>
*>          The factored form of the matrix A.  AFP contains the block<br>
*>          diagonal matrix D and the multipliers used to obtain the<br>
*>          factor U or L from the factorization A = U*D*U**T or<br>
*>          A = L*D*L**T as computed by SSPTRF, stored as a packed<br>
*>          triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by SSPTRF.<br>
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
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by SSPTRS.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> SSPSV computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspsv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspsv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspsv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPSV computes the solution to a real system of linear equations<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          See below for further details.<br>
*><br>
*>          On exit, the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by SSPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D, as<br>
*>          determined by SSPTRF.  If IPIV(k) > 0, then rows and columns<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realOTHERsolve<br>
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
	public void sspsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief <b> SSPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sspsvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sspsvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sspsvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X,<br>
*                          LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          FACT, UPLO<br>
*       INTEGER            INFO, LDB, LDX, N, NRHS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               AFP( * ), AP( * ), B( LDB, * ), BERR( * ),<br>
*      $                   FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPSVX uses the diagonal pivoting factorization A = U*D*U**T or<br>
*> A = L*D*L**T to compute the solution to a real system of linear<br>
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
*>          = 'F':  On entry, AFP and IPIV contain the factored form of<br>
*>                  A.  AP, AFP and IPIV will not be modified.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*>          AFP is REAL array, dimension<br>
*>                            (N*(N+1)/2)<br>
*>          If FACT = 'F', then AFP is an input argument and on entry<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by SSPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*><br>
*>          If FACT = 'N', then AFP is an output argument and on exit<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by SSPTRF, stored as<br>
*>          a packed triangular matrix in the same storage format as A.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          If FACT = 'F', then IPIV is an input argument and on entry<br>
*>          contains details of the interchanges and the block structure<br>
*>          of D, as determined by SSPTRF.<br>
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
*>          of D, as determined by SSPTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A.  If RCOND is less than the machine precision (in<br>
*>          particular, if RCOND = 0), the matrix is singular to working<br>
*>          precision.  This condition is indicated by a return code of<br>
*>          INFO > 0.<br>
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
*> \ingroup realOTHERsolve<br>
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
	public void sspsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SSPTRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPTRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPTRD( UPLO, N, AP, D, E, TAU, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               AP( * ), D( * ), E( * ), TAU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPTRD reduces a real symmetric matrix A stored in packed form to<br>
*> symmetric tridiagonal form T by an orthogonal similarity<br>
*> transformation: Q**T * A * Q = T.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the upper or lower triangle of the symmetric matrix<br>
*>          A, packed columnwise in a linear array.  The j-th column of A<br>
*>          is stored in the array AP as follows:<br>
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;<br>
*>          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.<br>
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal<br>
*>          of A are overwritten by the corresponding elements of the<br>
*>          tridiagonal matrix T, and the elements above the first<br>
*>          superdiagonal, with the array TAU, represent the orthogonal<br>
*>          matrix Q as a product of elementary reflectors; if UPLO<br>
*>          = 'L', the diagonal and first subdiagonal of A are over-<br>
*>          written by the corresponding elements of the tridiagonal<br>
*>          matrix T, and the elements below the first subdiagonal, with<br>
*>          the array TAU, represent the orthogonal matrix Q as a product<br>
*>          of elementary reflectors. See Further Details.<br>
*> \endverbatim<br>
*><br>
*> \param[out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is REAL array, dimension (N-1)<br>
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
*> \ingroup realOTHERcomputational<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
*>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,<br>
*>  overwriting A(i+2:n,i), and tau is stored in TAU(i).<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void ssptrd_(CHARACTER UPLO,INTEGER N,float[] AP,float[] D,float[] E,float[] TAU,INTEGER INFO);
/**
*> \brief \b SSPTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPTRF( UPLO, N, AP, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               AP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPTRF computes the factorization of a real symmetric matrix A stored<br>
*> in packed format using the Bunch-Kaufman diagonal pivoting method:<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
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
*> \ingroup realOTHERcomputational<br>
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
	public void ssptrf_(CHARACTER UPLO,INTEGER N,float[] AP,int[] IPIV,INTEGER INFO);
/**
*> \brief \b SSPTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPTRI( UPLO, N, AP, IPIV, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               AP( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPTRI computes the inverse of a real symmetric indefinite matrix<br>
*> A in packed storage using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by SSPTRF.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          On entry, the block diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by SSPTRF,<br>
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
*>          as determined by SSPTRF.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssptri_(CHARACTER UPLO,INTEGER N,float[] AP,int[] IPIV,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSPTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSPTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssptrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssptrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssptrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               AP( * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSPTRS solves a system of linear equations A*X = B with a real<br>
*> symmetric matrix A stored in packed format using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by SSPTRF.<br>
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
*>          AP is REAL array, dimension (N*(N+1)/2)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSPTRF, stored as a<br>
*>          packed triangular matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] IPIV<br>
*> \verbatim<br>
*>          IPIV is INTEGER array, dimension (N)<br>
*>          Details of the interchanges and the block structure of D<br>
*>          as determined by SSPTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b SSTEBZ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEBZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstebz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstebz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstebz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E,<br>
*                          M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,<br>
*                          INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          ORDER, RANGE<br>
*       INTEGER            IL, INFO, IU, M, N, NSPLIT<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )<br>
*       REAL               D( * ), E( * ), W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEBZ computes the eigenvalues of a symmetric tridiagonal<br>
*> matrix T.  The user may ask for all eigenvalues, all eigenvalues<br>
*> in the half-open interval (VL, VU], or the IL-th through IU-th<br>
*> eigenvalues.<br>
*><br>
*> To avoid overflow, the matrix must be scaled so that its<br>
*> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest<br>
*> accuracy, it should not be much smaller than that.<br>
*><br>
*> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal<br>
*> Matrix", Report CS41, Computer Science Dept., Stanford<br>
*> University, July 21, 1966.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] RANGE<br>
*> \verbatim<br>
*>          RANGE is CHARACTER*1<br>
*>          = 'A': ("All")   all eigenvalues will be found.<br>
*>          = 'V': ("Value") all eigenvalues in the half-open interval<br>
*>                           (VL, VU] will be found.<br>
*>          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the<br>
*>                           entire matrix) will be found.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ORDER<br>
*> \verbatim<br>
*>          ORDER is CHARACTER*1<br>
*>          = 'B': ("By Block") the eigenvalues will be grouped by<br>
*>                              split-off block (see IBLOCK, ISPLIT) and<br>
*>                              ordered from smallest to largest within<br>
*>                              the block.<br>
*>          = 'E': ("Entire matrix")<br>
*>                              the eigenvalues for the entire matrix<br>
*>                              will be ordered from smallest to<br>
*>                              largest.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the tridiagonal matrix T.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues.  Eigenvalues less than or equal<br>
*>          to VL, or greater than VU, will not be returned.  VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
*><br>
*>          If RANGE='V', the upper bound of the interval to<br>
*>          be searched for eigenvalues.  Eigenvalues less than or equal<br>
*>          to VL, or greater than VU, will not be returned.  VL < VU.<br>
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
*>          ABSTOL is REAL<br>
*>          The absolute tolerance for the eigenvalues.  An eigenvalue<br>
*>          (or cluster) is considered to be located if it has been<br>
*>          determined to lie in an interval whose width is ABSTOL or<br>
*>          less.  If ABSTOL is less than or equal to zero, then ULP*|T|<br>
*>          will be used, where |T| means the 1-norm of T.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The (n-1) off-diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The actual number of eigenvalues found. 0 <= M <= N.<br>
*>          (See also the description of INFO=2,3.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] NSPLIT<br>
*> \verbatim<br>
*>          NSPLIT is INTEGER<br>
*>          The number of diagonal blocks in the matrix T.<br>
*>          1 <= NSPLIT <= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] W<br>
*> \verbatim<br>
*>          W is REAL array, dimension (N)<br>
*>          On exit, the first M elements of W will contain the<br>
*>          eigenvalues.  (SSTEBZ may use the remaining N-M elements as<br>
*>          workspace.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IBLOCK<br>
*> \verbatim<br>
*>          IBLOCK is INTEGER array, dimension (N)<br>
*>          At each row/column j where E(j) is zero or small, the<br>
*>          matrix T is considered to split into a block diagonal<br>
*>          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which<br>
*>          block (from 1 to the number of blocks) the eigenvalue W(i)<br>
*>          belongs.  (SSTEBZ may use the remaining N-M elements as<br>
*>          workspace.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into submatrices.<br>
*>          The first submatrix consists of rows/columns 1 to ISPLIT(1),<br>
*>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),<br>
*>          etc., and the NSPLIT-th consists of rows/columns<br>
*>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.<br>
*>          (Only the first NSPLIT elements will actually be used, but<br>
*>          since the user cannot know a priori what value NSPLIT will<br>
*>          have, N words must be reserved for ISPLIT.)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (4*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (3*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  some or all of the eigenvalues failed to converge or<br>
*>                were not computed:<br>
*>                =1 or 3: Bisection failed to converge for some<br>
*>                        eigenvalues; these eigenvalues are flagged by a<br>
*>                        negative block number.  The effect is that the<br>
*>                        eigenvalues may not be as accurate as the<br>
*>                        absolute and relative tolerances.  This is<br>
*>                        generally caused by unexpectedly inaccurate<br>
*>                        arithmetic.<br>
*>                =2 or 3: RANGE='I' only: Not all of the eigenvalues<br>
*>                        IL:IU were found.<br>
*>                        Effect: M < IU+1-IL<br>
*>                        Cause:  non-monotonic arithmetic, causing the<br>
*>                                Sturm sequence to be non-monotonic.<br>
*>                        Cure:   recalculate, using RANGE='A', and pick<br>
*>                                out eigenvalues IL:IU.  In some cases,<br>
*>                                increasing the PARAMETER "FUDGE" may<br>
*>                                make things work.<br>
*>                = 4:    RANGE='I', and the Gershgorin interval<br>
*>                        initially used was too small.  No eigenvalues<br>
*>                        were computed.<br>
*>                        Probable cause: your machine has sloppy<br>
*>                                        floating-point arithmetic.<br>
*>                        Cure: Increase the PARAMETER "FUDGE",<br>
*>                              recompile, and try again.<br>
*> \endverbatim<br>
*<br>
*> \par Internal Parameters:<br>
*  =========================<br>
*><br>
*> \verbatim<br>
*>  RELFAC  REAL, default = 2.0e0<br>
*>          The relative tolerance.  An interval (a,b] lies within<br>
*>          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),<br>
*>          where "ulp" is the machine precision (distance from 1 to<br>
*>          the next larger floating point number.)<br>
*><br>
*>  FUDGE   REAL, default = 2<br>
*>          A "fudge factor" to widen the Gershgorin intervals.  Ideally,<br>
*>          a value of 1 should work, but on machines with sloppy<br>
*>          arithmetic, this needs to be larger.  The default for<br>
*>          publicly released versions should be large enough to handle<br>
*>          the worst machine around.  Note that this has no effect<br>
*>          on accuracy of the solution.<br>
*> \endverbatim<br>
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
*  =====================================================================<br>
*/
	public void sstebz_(CHARACTER RANGE,CHARACTER ORDER,INTEGER N,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,float[] D,float[] E,INTEGER M,INTEGER NSPLIT,float[] W,int[] IBLOCK,int[] ISPLIT,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SSTEDC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEDC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstedc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstedc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstedc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,<br>
*                          LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPZ<br>
*       INTEGER            INFO, LDZ, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEDC computes all eigenvalues and, optionally, eigenvectors of a<br>
*> symmetric tridiagonal matrix using the divide and conquer method.<br>
*> The eigenvectors of a full or band real symmetric matrix can also be<br>
*> found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this<br>
*> matrix to tridiagonal form.<br>
*><br>
*> This code makes very mild assumptions about floating point<br>
*> arithmetic. It will work on machines with a guard digit in<br>
*> add/subtract, or on those binary machines without guard digits<br>
*> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.<br>
*> It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.  See SLAED3 for details.<br>
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
*>          = 'V':  Compute eigenvectors of original dense symmetric<br>
*>                  matrix also.  On entry, Z contains the orthogonal<br>
*>                  matrix used to reduce the original matrix to<br>
*>                  tridiagonal form.<br>
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
*>          D is REAL array, dimension (N)<br>
*>          On entry, the diagonal elements of the tridiagonal matrix.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          On entry, the subdiagonal elements of the tridiagonal matrix.<br>
*>          On exit, E has been destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ,N)<br>
*>          On entry, if COMPZ = 'V', then Z contains the orthogonal<br>
*>          matrix used in the reduction to tridiagonal form.<br>
*>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the<br>
*>          orthonormal eigenvectors of the original symmetric matrix,<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.<br>
*>          If COMPZ = 'V' and N > 1 then LWORK must be at least<br>
*>                         ( 1 + 3*N + 2*N*lg N + 4*N**2 ),<br>
*>                         where lg( N ) = smallest integer k such<br>
*>                         that 2**k >= N.<br>
*>          If COMPZ = 'I' and N > 1 then LWORK must be at least<br>
*>                         ( 1 + 4*N + N**2 ).<br>
*>          Note that for COMPZ = 'I' or 'V', then if N is less than or<br>
*>          equal to the minimum divide size, usually 25, then LWORK need<br>
*>          only be max(1,2*(N-1)).<br>
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
*>          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.<br>
*>          If COMPZ = 'V' and N > 1 then LIWORK must be at least<br>
*>                         ( 6 + 6*N + 5*N*lg N ).<br>
*>          If COMPZ = 'I' and N > 1 then LIWORK must be at least<br>
*>                         ( 3 + 5*N ).<br>
*>          Note that for COMPZ = 'I' or 'V', then if N is less than or<br>
*>          equal to the minimum divide size, usually 25, then LIWORK<br>
*>          need only be 1.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA \n<br>
*>  Modified by Francoise Tisseur, University of Tennessee<br>
*><br>
*  =====================================================================<br>
*/
	public void sstedc_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b SSTEGR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEGR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstegr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstegr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstegr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,<br>
*                  ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,<br>
*                  LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE<br>
*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N<br>
*       REAL             ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       REAL               D( * ), E( * ), W( * ), WORK( * )<br>
*       REAL               Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEGR computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has<br>
*> a well defined set of pairwise different real eigenvalues, the corresponding<br>
*> real eigenvectors are pairwise orthogonal.<br>
*><br>
*> The spectrum may be computed either completely or partially by specifying<br>
*> either an interval (VL,VU] or a range of indices IL:IU for the desired<br>
*> eigenvalues.<br>
*><br>
*> SSTEGR is a compatibility wrapper around the improved SSTEMR routine.<br>
*> See SSTEMR for further details.<br>
*><br>
*> One important change is that the ABSTOL parameter no longer provides any<br>
*> benefit and hence is no longer used.<br>
*><br>
*> Note : SSTEGR and SSTEMR work only on machines which follow<br>
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
*>          D is REAL array, dimension (N)<br>
*>          On entry, the N diagonal elements of the tridiagonal matrix<br>
*>          T. On exit, D is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          On entry, the (N-1) subdiagonal elements of the tridiagonal<br>
*>          matrix T in elements 1 to N-1 of E. E(N) need not be set on<br>
*>          input, but is used internally as workspace.<br>
*>          On exit, E is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M) )<br>
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
*>          WORK is REAL array, dimension (LWORK)<br>
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
*>          > 0:  if INFO = 1X, internal error in SLARRE,<br>
*>                if INFO = 2X, internal error in SLARRV.<br>
*>                Here, the digit X = ABS( IINFO ) < 10, where IINFO is<br>
*>                the nonzero error code returned by SLARRE or<br>
*>                SLARRV, respectively.<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Inderjit Dhillon, IBM Almaden, USA \n<br>
*> Osni Marques, LBNL/NERSC, USA \n<br>
*> Christof Voemel, LBNL/NERSC, USA \n<br>
*<br>
*  =====================================================================<br>
*/
	public void sstegr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b SSTEIN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEIN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstein.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstein.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstein.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,<br>
*                          IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDZ, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),<br>
*      $                   IWORK( * )<br>
*       REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEIN computes the eigenvectors of a real symmetric tridiagonal<br>
*> matrix T corresponding to specified eigenvalues, using inverse<br>
*> iteration.<br>
*><br>
*> The maximum number of iterations allowed for each eigenvector is<br>
*> specified by an internal parameter MAXITS (currently set to 5).<br>
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
*>          D is REAL array, dimension (N)<br>
*>          The n diagonal elements of the tridiagonal matrix T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The (n-1) subdiagonal elements of the tridiagonal matrix<br>
*>          T, in elements 1 to N-1.<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements of W contain the eigenvalues for<br>
*>          which eigenvectors are to be computed.  The eigenvalues<br>
*>          should be grouped by split-off block and ordered from<br>
*>          smallest to largest within the block.  ( The output array<br>
*>          W from SSTEBZ with ORDER = 'B' is expected here. )<br>
*> \endverbatim<br>
*><br>
*> \param[in] IBLOCK<br>
*> \verbatim<br>
*>          IBLOCK is INTEGER array, dimension (N)<br>
*>          The submatrix indices associated with the corresponding<br>
*>          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to<br>
*>          the first submatrix from the top, =2 if W(i) belongs to<br>
*>          the second submatrix, etc.  ( The output array IBLOCK<br>
*>          from SSTEBZ is expected here. )<br>
*> \endverbatim<br>
*><br>
*> \param[in] ISPLIT<br>
*> \verbatim<br>
*>          ISPLIT is INTEGER array, dimension (N)<br>
*>          The splitting points, at which T breaks up into submatrices.<br>
*>          The first submatrix consists of rows/columns 1 to<br>
*>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1<br>
*>          through ISPLIT( 2 ), etc.<br>
*>          ( The output array ISPLIT from SSTEBZ is expected here. )<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, M)<br>
*>          The computed eigenvectors.  The eigenvector associated<br>
*>          with the eigenvalue W(i) is stored in the i-th column of<br>
*>          Z.  Any vector which fails to converge is set to its current<br>
*>          iterate after MAXITS iterations.<br>
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
*>          WORK is REAL array, dimension (5*N)<br>
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
*>          = 0: successful exit.<br>
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
*> \ingroup realOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sstein_(INTEGER N,float[] D,float[] E,INTEGER M,float[] W,int[] IBLOCK,int[] ISPLIT,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSTEMR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEMR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstemr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstemr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstemr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,<br>
*                          M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK,<br>
*                          IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE<br>
*       LOGICAL            TRYRAC<br>
*       INTEGER            IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N<br>
*       REAL               VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       REAL               D( * ), E( * ), W( * ), WORK( * )<br>
*       REAL               Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEMR computes selected eigenvalues and, optionally, eigenvectors<br>
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
*> 1.SSTEMR works only on machines which follow IEEE-754<br>
*> floating-point standard in their handling of infinities and NaNs.<br>
*> This permits the use of efficient inner loops avoiding a check for<br>
*> zero divisors.<br>
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
*>          D is REAL array, dimension (N)<br>
*>          On entry, the N diagonal elements of the tridiagonal matrix<br>
*>          T. On exit, D is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
*>          On entry, the (N-1) subdiagonal elements of the tridiagonal<br>
*>          matrix T in elements 1 to N-1 of E. E(N) need not be set on<br>
*>          input, but is used internally as workspace.<br>
*>          On exit, E is overwritten.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*><br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M) )<br>
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
*>          WORK is REAL array, dimension (LWORK)<br>
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
*>          > 0:  if INFO = 1X, internal error in SLARRE,<br>
*>                if INFO = 2X, internal error in SLARRV.<br>
*>                Here, the digit X = ABS( IINFO ) < 10, where IINFO is<br>
*>                the nonzero error code returned by SLARRE or<br>
*>                SLARRV, respectively.<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Beresford Parlett, University of California, Berkeley, USA \n<br>
*> Jim Demmel, University of California, Berkeley, USA \n<br>
*> Inderjit Dhillon, University of Texas, Austin, USA \n<br>
*> Osni Marques, LBNL/NERSC, USA \n<br>
*> Christof Voemel, University of California, Berkeley, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void sstemr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,INTEGER M,float[] W,float[] Z,INTEGER LDZ,INTEGER NZC,int[] ISUPPZ,LOGICAL TRYRAC,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b SSTEQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssteqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssteqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssteqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPZ<br>
*       INTEGER            INFO, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEQR computes all eigenvalues and, optionally, eigenvectors of a<br>
*> symmetric tridiagonal matrix using the implicit QL or QR method.<br>
*> The eigenvectors of a full or band symmetric matrix can also be found<br>
*> if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to<br>
*> tridiagonal form.<br>
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
*>                  symmetric matrix.  On entry, Z must contain the<br>
*>                  orthogonal matrix used to reduce the original matrix<br>
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
*>          D is REAL array, dimension (N)<br>
*>          On entry, the diagonal elements of the tridiagonal matrix.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
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
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          On entry, if  COMPZ = 'V', then Z contains the orthogonal<br>
*>          matrix used in the reduction to tridiagonal form.<br>
*>          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the<br>
*>          orthonormal eigenvectors of the original symmetric matrix,<br>
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
*>          WORK is REAL array, dimension (max(1,2*N-2))<br>
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
*>                matrix which is orthogonally similar to the original<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssteqr_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSTERF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTERF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssterf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssterf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssterf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTERF( N, D, E, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTERF computes all eigenvalues of a symmetric tridiagonal matrix<br>
*> using the Pal-Walker-Kahan variant of the QL or QR algorithm.<br>
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
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
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
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  the algorithm failed to find all of the eigenvalues in<br>
*>                a total of 30*N iterations; if INFO = i, then i<br>
*>                elements of E have not converged to zero.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
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
	public void ssterf_(INTEGER N,float[] D,float[] E,INTEGER INFO);
/**
*> \brief <b> SSTEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ<br>
*       INTEGER            INFO, LDZ, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEV computes all eigenvalues and, optionally, eigenvectors of a<br>
*> real symmetric tridiagonal matrix A.<br>
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
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix<br>
*>          A.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix A, stored in elements 1 to N-1 of E.<br>
*>          On exit, the contents of E are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal<br>
*>          eigenvectors of the matrix A, with the i-th column of Z<br>
*>          holding the eigenvector associated with D(i).<br>
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
*>          WORK is REAL array, dimension (max(1,2*N-2))<br>
*>          If JOBZ = 'N', WORK is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of E did not converge to zero.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2011<br>
*<br>
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sstev_(CHARACTER JOBZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
/**
*> \brief <b> SSTEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,<br>
*                          LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ<br>
*       INTEGER            INFO, LDZ, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEVD computes all eigenvalues and, optionally, eigenvectors of a<br>
*> real symmetric tridiagonal matrix. If eigenvectors are desired, it<br>
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
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix<br>
*>          A.<br>
*>          On exit, if INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix A, stored in elements 1 to N-1 of E.<br>
*>          On exit, the contents of E are destroyed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, N)<br>
*>          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal<br>
*>          eigenvectors of the matrix A, with the i-th column of Z<br>
*>          holding the eigenvector associated with D(i).<br>
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
*>          WORK is REAL array,<br>
*>                                         dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.<br>
*>          If JOBZ  = 'V' and N > 1 then LWORK must be at least<br>
*>                         ( 1 + 4*N + N**2 ).<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
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
*>          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.<br>
*>          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, the algorithm failed to converge; i<br>
*>                off-diagonal elements of E did not converge to zero.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date November 2011<br>
*<br>
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sstevd_(CHARACTER JOBZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> SSTEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEVR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,<br>
*                          M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,<br>
*                          LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE<br>
*       INTEGER            IL, INFO, IU, LDZ, LIWORK, LWORK, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEVR computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric tridiagonal matrix T.  Eigenvalues and<br>
*> eigenvectors can be selected by specifying either a range of values<br>
*> or a range of indices for the desired eigenvalues.<br>
*><br>
*> Whenever possible, SSTEVR calls SSTEMR to compute the<br>
*> eigenspectrum using Relatively Robust Representations.  SSTEMR<br>
*> computes eigenvalues by the dqds algorithm, while orthogonal<br>
*> eigenvectors are computed from various "good" L D L^T representations<br>
*> (also known as Relatively Robust Representations). Gram-Schmidt<br>
*> orthogonalization is avoided as far as possible. More specifically,<br>
*> the various steps of the algorithm are as follows. For the i-th<br>
*> unreduced block of T,<br>
*>    (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T<br>
*>         is a relatively robust representation,<br>
*>    (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high<br>
*>        relative accuracy by the dqds algorithm,<br>
*>    (c) If there is a cluster of close eigenvalues, "choose" sigma_i<br>
*>        close to the cluster, and go to step (a),<br>
*>    (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,<br>
*>        compute the corresponding eigenvector by forming a<br>
*>        rank-revealing twisted factorization.<br>
*> The desired accuracy of the output can be specified by the input<br>
*> parameter ABSTOL.<br>
*><br>
*> For more details, see "A new O(n^2) algorithm for the symmetric<br>
*> tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,<br>
*> Computer Science Division Technical Report No. UCB//CSD-97-971,<br>
*> UC Berkeley, May 1997.<br>
*><br>
*><br>
*> Note 1 : SSTEVR calls SSTEMR when the full spectrum is requested<br>
*> on machines which conform to the ieee-754 floating point standard.<br>
*> SSTEVR calls SSTEBZ and SSTEIN on non-ieee machines and<br>
*> when partial spectrum requests are made.<br>
*><br>
*> Normal execution of SSTEMR may create NaNs and infinities and<br>
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
*>          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and<br>
*>          SSTEIN are called<br>
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
*>          On entry, the n diagonal elements of the tridiagonal matrix<br>
*>          A.<br>
*>          On exit, D may be multiplied by a constant factor chosen<br>
*>          to avoid over/underflow in computing the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (max(1,N-1))<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix A in elements 1 to N-1 of E.<br>
*>          On exit, E may be multiplied by a constant factor chosen<br>
*>          to avoid over/underflow in computing the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that<br>
*>          eigenvalues are computed to high relative accuracy when<br>
*>          possible in future releases.  The current code does not<br>
*>          make any guarantees about high relative accuracy, but<br>
*>          future releases will. See J. Barlow and J. Demmel,<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M) )<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal (and<br>
*>          minimal) LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.  LWORK >= 20*N.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if INFO = 0, IWORK(1) returns the optimal (and<br>
*>          minimal) LIWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.  LIWORK >= 10*N.<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
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
*> \ingroup realOTHEReigen<br>
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
	public void sstevr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> SSTEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSTEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,<br>
*                          M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE<br>
*       INTEGER            IL, INFO, IU, LDZ, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSTEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric tridiagonal matrix A.  Eigenvalues and<br>
*> eigenvectors can be selected by specifying either a range of values<br>
*> or a range of indices for the desired eigenvalues.<br>
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
*>          D is REAL array, dimension (N)<br>
*>          On entry, the n diagonal elements of the tridiagonal matrix<br>
*>          A.<br>
*>          On exit, D may be multiplied by a constant factor chosen<br>
*>          to avoid over/underflow in computing the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (max(1,N-1))<br>
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal<br>
*>          matrix A in elements 1 to N-1 of E.<br>
*>          On exit, E may be multiplied by a constant factor chosen<br>
*>          to avoid over/underflow in computing the eigenvalues.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VL<br>
*> \verbatim<br>
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
*>          The absolute error tolerance for the eigenvalues.<br>
*>          An approximate eigenvalue is accepted as converged<br>
*>          when it is determined to lie in an interval [a,b]<br>
*>          of width less than or equal to<br>
*><br>
*>                  ABSTOL + EPS *   max( |a|,|b| ) ,<br>
*><br>
*>          where EPS is the machine precision.  If ABSTOL is less<br>
*>          than or equal to zero, then  EPS*|T|  will be used in<br>
*>          its place, where |T| is the 1-norm of the tridiagonal<br>
*>          matrix.<br>
*><br>
*>          Eigenvalues will be computed most accurately when ABSTOL is<br>
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M) )<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
*>          corresponding to the selected eigenvalues, with the i-th<br>
*>          column of Z holding the eigenvector associated with W(i).<br>
*>          If an eigenvector fails to converge (INFO > 0), then that<br>
*>          column of Z contains the latest approximation to the<br>
*>          eigenvector, and the index of the eigenvector is returned<br>
*>          in IFAIL.  If JOBZ = 'N', then Z is not referenced.<br>
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
*>          WORK is REAL array, dimension (5*N)<br>
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
*> \ingroup realOTHEReigen<br>
*<br>
*  =====================================================================<br>
*/
	public void sstevx_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSYCON<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYCON + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssycon.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssycon.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssycon.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYCON estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a real symmetric matrix A using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by SSYTRF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
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
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssycon_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SSYCONV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYCONV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyconv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyconv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyconv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, WAY<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), E( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYCONV convert A given by TRF into L and D and vice-versa.<br>
*> Get Non-diag elements of D (returned in workspace) and <br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssyconv_(CHARACTER UPLO,CHARACTER WAY,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] E,INTEGER INFO);
/**
*> \brief \b SSYCON_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYCON_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssycon_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssycon_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssycon_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYCON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND,<br>
*                               WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               ANORM, RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYCON_ROOK estimates the reciprocal of the condition number (in the<br>
*> 1-norm) of a real symmetric matrix A using the factorization<br>
*> A = U*D*U**T or A = L*D*L**T computed by SSYTRF_ROOK.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSYTRF_ROOK.<br>
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
*>          as determined by SSYTRF_ROOK.<br>
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
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an<br>
*>          estimate of the 1-norm of inv(A) computed in this routine.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (2*N)<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup realSYcomputational<br>
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
	public void ssycon_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SSYEQUB<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYEQUB + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyequb.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyequb.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyequb.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INFO, LDA, N<br>
*       REAL               AMAX, SCOND<br>
*       CHARACTER          UPLO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), S( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYEQUB computes row and column scalings intended to equilibrate a<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>  Livne, O.E. and Golub, G.H., "Scaling by Binormalization", \n<br>
*>  Numerical Algorithms, vol. 35, no. 1, pp. 97-120, January 2004. \n<br>
*>  DOI 10.1023/B:NUMA.0000016606.32820.69 \n<br>
*>  Tech report version: http://ruready.utah.edu/archive/papers/bin.pdf<br>
*><br>
*  =====================================================================<br>
*/
	public void ssyequb_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,float[] WORK,INTEGER INFO);
/**
*> \brief <b> SSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYEV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyev.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyev.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyev.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYEV computes all eigenvalues and, optionally, eigenvectors of a<br>
*> real symmetric matrix A.<br>
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
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          The length of the array WORK.  LWORK >= max(1,3*N-1).<br>
*>          For optimal efficiency, LWORK >= (NB+2)*N,<br>
*>          where NB is the blocksize for SSYTRD returned by ILAENV.<br>
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
*> \ingroup realSYeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssyev_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] W,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> SSYEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYEVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK,<br>
*                          LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, LDA, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYEVD computes all eigenvalues and, optionally, eigenvectors of a<br>
*> real symmetric matrix A. If eigenvectors are desired, it uses a<br>
*> divide and conquer algorithm.<br>
*><br>
*> The divide and conquer algorithm makes very mild assumptions about<br>
*> floating point arithmetic. It will work on machines with a guard<br>
*> digit in add/subtract, or on those binary machines without guard<br>
*> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or<br>
*> Cray-2. It could conceivably fail on hexadecimal or decimal machines<br>
*> without guard digits, but we know of none.<br>
*><br>
*> Because of large use of BLAS of level 3, SSYEVD needs N**2 more<br>
*> workspace than SSYEVX.<br>
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
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array,<br>
*>                                         dimension (LWORK)<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.<br>
*>          If N <= 1,               LWORK must be at least 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.<br>
*>          If JOBZ = 'V' and N > 1, LWORK must be at least <br>
*>                                                1 + 6*N + 2*N**2.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
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
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
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
*> \ingroup realSYeigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Jeff Rutter, Computer Science Division, University of California<br>
*> at Berkeley, USA \n<br>
*>  Modified by Francoise Tisseur, University of Tennessee \n<br>
*>  Modified description of INFO. Sven, 16 Feb 05. \n<br>
*><br>
*  =====================================================================<br>
*/
	public void ssyevd_(CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] W,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> SSYEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYEVR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,<br>
*                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,<br>
*                          IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            ISUPPZ( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYEVR computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric matrix A.  Eigenvalues and eigenvectors can be<br>
*> selected by specifying either a range of values or a range of<br>
*> indices for the desired eigenvalues.<br>
*><br>
*> SSYEVR first reduces the matrix A to tridiagonal form T with a call<br>
*> to SSYTRD.  Then, whenever possible, SSYEVR calls SSTEMR to compute<br>
*> the eigenspectrum using Relatively Robust Representations.  SSTEMR<br>
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
*> For more details, see SSTEMR's documentation and:<br>
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
*> Note 1 : SSYEVR calls SSTEMR when the full spectrum is requested<br>
*> on machines which conform to the ieee-754 floating point standard.<br>
*> SSYEVR calls SSTEBZ and SSTEIN on non-ieee machines and<br>
*> when partial spectrum requests are made.<br>
*><br>
*> Normal execution of SSTEMR may create NaNs and infinities and<br>
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
*>          For RANGE = 'V' or 'I' and IU - IL < N - 1, SSTEBZ and<br>
*>          SSTEIN are called<br>
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
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
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
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          SLAMCH( 'Safe minimum' ).  Doing so will guarantee that<br>
*>          eigenvalues are computed to high relative accuracy when<br>
*>          possible in future releases.  The current code does not<br>
*>          make any guarantees about high relative accuracy, but<br>
*>          future releases will. See J. Barlow and J. Demmel,<br>
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
*>          W is REAL array, dimension (N)<br>
*>          The first M elements contain the selected eigenvalues in<br>
*>          ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M))<br>
*>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z<br>
*>          contain the orthonormal eigenvectors of the matrix A<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.  LWORK >= max(1,26*N).<br>
*>          For optimal efficiency, LWORK >= (NB+6)*N,<br>
*>          where NB is the max of the blocksize for SSYTRD and SORMTR<br>
*>          returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IWORK<br>
*> \verbatim<br>
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))<br>
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LIWORK<br>
*> \verbatim<br>
*>          LIWORK is INTEGER<br>
*>          The dimension of the array IWORK.  LIWORK >= max(1,10*N).<br>
*><br>
*>          If LIWORK = -1, then a workspace query is assumed; the<br>
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
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
*> \ingroup realSYeigen<br>
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
	public void ssyevr_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief <b> SSYEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYEVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyevx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyevx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyevx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,<br>
*                          ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,<br>
*                          IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYEVX computes selected eigenvalues and, optionally, eigenvectors<br>
*> of a real symmetric matrix A.  Eigenvalues and eigenvectors can be<br>
*> selected by specifying either a range of values or a range of indices<br>
*> for the desired eigenvalues.<br>
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
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
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
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.<br>
*>          If this routine returns with INFO>0, indicating that some<br>
*>          eigenvectors did not converge, try setting ABSTOL to<br>
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          On normal exit, the first M elements contain the selected<br>
*>          eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M))<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of the array WORK.  LWORK >= 1, when N <= 1;<br>
*>          otherwise 8*N.<br>
*>          For optimal efficiency, LWORK >= (NB+3)*N,<br>
*>          where NB is the max of the blocksize for SSYTRD and SORMTR<br>
*>          returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*> \ingroup realSYeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssyevx_(CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSYGS2 reduces a symmetric definite generalized eigenproblem to standard form, using the factorization results obtained from spotrf (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYGS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygs2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygs2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygs2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, N<br>
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
*> SSYGS2 reduces a real symmetric-definite generalized eigenproblem<br>
*> to standard form.<br>
*><br>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,<br>
*> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)<br>
*><br>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or<br>
*> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T *A*L.<br>
*><br>
*> B must have been previously factorized as U**T *U or L*L**T by SPOTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);<br>
*>          = 2 or 3: compute U*A*U**T or L**T *A*L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          Specifies whether the upper or lower triangular part of the<br>
*>          symmetric matrix A is stored, and how B has been factorized.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
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
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,N)<br>
*>          The triangular factor from the Cholesky factorization of B,<br>
*>          as returned by SPOTRF.<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssygs2_(INTEGER ITYPE,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b SSYGST<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYGST + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygst.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygst.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygst.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, N<br>
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
*> SSYGST reduces a real symmetric-definite generalized eigenproblem<br>
*> to standard form.<br>
*><br>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,<br>
*> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)<br>
*><br>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or<br>
*> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.<br>
*><br>
*> B must have been previously factorized as U**T*U or L*L**T by SPOTRF.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ITYPE<br>
*> \verbatim<br>
*>          ITYPE is INTEGER<br>
*>          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);<br>
*>          = 2 or 3: compute U*A*U**T or L**T*A*L.<br>
*> \endverbatim<br>
*><br>
*> \param[in] UPLO<br>
*> \verbatim<br>
*>          UPLO is CHARACTER*1<br>
*>          = 'U':  Upper triangle of A is stored and B is factored as<br>
*>                  U**T*U;<br>
*>          = 'L':  Lower triangle of A is stored and B is factored as<br>
*>                  L*L**T.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
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
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,N)<br>
*>          The triangular factor from the Cholesky factorization of B,<br>
*>          as returned by SPOTRF.<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssygst_(INTEGER ITYPE,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b SSYGV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYGV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,<br>
*                         LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYGV computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a real generalized symmetric-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.<br>
*> Here A and B are assumed to be symmetric and B is also<br>
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
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*><br>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the<br>
*>          matrix Z of eigenvectors.  The eigenvectors are normalized<br>
*>          as follows:<br>
*>          if ITYPE = 1 or 2, Z**T*B*Z = I;<br>
*>          if ITYPE = 3, Z**T*inv(B)*Z = I.<br>
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
*>          B is REAL array, dimension (LDB, N)<br>
*>          On entry, the symmetric positive definite matrix B.<br>
*>          If UPLO = 'U', the leading N-by-N upper triangular part of B<br>
*>          contains the upper triangular part of the matrix B.<br>
*>          If UPLO = 'L', the leading N-by-N lower triangular part of B<br>
*>          contains the lower triangular part of the matrix B.<br>
*><br>
*>          On exit, if INFO <= N, the part of B containing the matrix is<br>
*>          overwritten by the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**T*U or B = L*L**T.<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          The length of the array WORK.  LWORK >= max(1,3*N-1).<br>
*>          For optimal efficiency, LWORK >= (NB+2)*N,<br>
*>          where NB is the blocksize for SSYTRD returned by ILAENV.<br>
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
*>          > 0:  SPOTRF or SSYEV returned an error code:<br>
*>             <= N:  if INFO = i, SSYEV failed to converge;<br>
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
*> \ingroup realSYeigen<br>
*<br>
*  =====================================================================<br>
*/
	public void ssygv_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] W,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SSYGVD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYGVD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygvd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygvd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygvd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,<br>
*                          LWORK, IWORK, LIWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, UPLO<br>
*       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IWORK( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYGVD computes all the eigenvalues, and optionally, the eigenvectors<br>
*> of a real generalized symmetric-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and<br>
*> B are assumed to be symmetric and B is also positive definite.<br>
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
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*><br>
*>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the<br>
*>          matrix Z of eigenvectors.  The eigenvectors are normalized<br>
*>          as follows:<br>
*>          if ITYPE = 1 or 2, Z**T*B*Z = I;<br>
*>          if ITYPE = 3, Z**T*inv(B)*Z = I.<br>
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
*>          B is REAL array, dimension (LDB, N)<br>
*>          On entry, the symmetric matrix B.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of B contains the<br>
*>          upper triangular part of the matrix B.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of B contains<br>
*>          the lower triangular part of the matrix B.<br>
*><br>
*>          On exit, if INFO <= N, the part of B containing the matrix is<br>
*>          overwritten by the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**T*U or B = L*L**T.<br>
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
*>          W is REAL array, dimension (N)<br>
*>          If INFO = 0, the eigenvalues in ascending order.<br>
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
*>          If N <= 1,               LWORK >= 1.<br>
*>          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.<br>
*>          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal sizes of the WORK and IWORK<br>
*>          arrays, returns these values as the first entries of the WORK<br>
*>          and IWORK arrays, and no error message related to LWORK or<br>
*>          LIWORK is issued by XERBLA.<br>
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
*>          routine only calculates the optimal sizes of the WORK and<br>
*>          IWORK arrays, returns these values as the first entries of<br>
*>          the WORK and IWORK arrays, and no error message related to<br>
*>          LWORK or LIWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  SPOTRF or SSYEVD returned an error code:<br>
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
*> \ingroup realSYeigen<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Modified so that no backsubstitution is performed if SSYEVD fails to<br>
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
	public void ssygvd_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] W,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
/**
*> \brief \b SSYGVX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYGVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,<br>
*                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,<br>
*                          LWORK, IWORK, IFAIL, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOBZ, RANGE, UPLO<br>
*       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N<br>
*       REAL               ABSTOL, VL, VU<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IFAIL( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYGVX computes selected eigenvalues, and optionally, eigenvectors<br>
*> of a real generalized symmetric-definite eigenproblem, of the form<br>
*> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A<br>
*> and B are assumed to be symmetric and B is also positive definite.<br>
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
*>          = 'U':  Upper triangle of A and B are stored;<br>
*>          = 'L':  Lower triangle of A and B are stored.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix pencil (A,B).  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of A contains the<br>
*>          upper triangular part of the matrix A.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of A contains<br>
*>          the lower triangular part of the matrix A.<br>
*><br>
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
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDA, N)<br>
*>          On entry, the symmetric matrix B.  If UPLO = 'U', the<br>
*>          leading N-by-N upper triangular part of B contains the<br>
*>          upper triangular part of the matrix B.  If UPLO = 'L',<br>
*>          the leading N-by-N lower triangular part of B contains<br>
*>          the lower triangular part of the matrix B.<br>
*><br>
*>          On exit, if INFO <= N, the part of B containing the matrix is<br>
*>          overwritten by the triangular factor U or L from the Cholesky<br>
*>          factorization B = U**T*U or B = L*L**T.<br>
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
*>          VL is REAL<br>
*>          If RANGE='V', the lower bound of the interval to<br>
*>          be searched for eigenvalues. VL < VU.<br>
*>          Not referenced if RANGE = 'A' or 'I'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] VU<br>
*> \verbatim<br>
*>          VU is REAL<br>
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
*>          ABSTOL is REAL<br>
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
*>          2*SLAMCH('S').<br>
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
*>          W is REAL array, dimension (N)<br>
*>          On normal exit, the first M elements contain the selected<br>
*>          eigenvalues in ascending order.<br>
*> \endverbatim<br>
*><br>
*> \param[out] Z<br>
*> \verbatim<br>
*>          Z is REAL array, dimension (LDZ, max(1,M))<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of the array WORK.  LWORK >= max(1,8*N).<br>
*>          For optimal efficiency, LWORK >= (NB+3)*N,<br>
*>          where NB is the blocksize for SSYTRD returned by ILAENV.<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
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
*>          > 0:  SPOTRF or SSYEVX returned an error code:<br>
*>             <= N:  if INFO = i, SSYEVX failed to converge;<br>
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
*> \ingroup realSYeigen<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA<br>
*<br>
*  =====================================================================<br>
*/
	public void ssygvx_(INTEGER ITYPE,CHARACTER JOBZ,CHARACTER RANGE,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
/**
*> \brief \b SSYRFS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYRFS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyrfs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyrfs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyrfs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,<br>
*                          X, LDX, FERR, BERR, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYRFS improves the computed solution to a system of linear<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>          The factored form of the matrix A.  AF contains the block<br>
*>          diagonal matrix D and the multipliers used to obtain the<br>
*>          factor U or L from the factorization A = U*D*U**T or<br>
*>          A = L*D*L**T as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
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
*> \param[in,out] X<br>
*> \verbatim<br>
*>          X is REAL array, dimension (LDX,NRHS)<br>
*>          On entry, the solution matrix X, as computed by SSYTRS.<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssyrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SSYRFSX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYRFSX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyrfsx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyrfsx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyrfsx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYRFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                           S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS,<br>
*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,<br>
*                           WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO, EQUED<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX, * ), WORK( * )<br>
*       REAL               S( * ), PARAMS( * ), BERR( * ),<br>
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
*>    SSYRFSX improves the computed solution to a system of linear<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     The factored form of the matrix A.  AF contains the block<br>
*>     diagonal matrix D and the multipliers used to obtain the<br>
*>     factor U or L from the factorization A = U*D*U**T or A =<br>
*>     L*D*L**T as computed by SSYTRF.<br>
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
*>     as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] S<br>
*> \verbatim<br>
*>          S is REAL array, dimension (N)<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          WORK is REAL array, dimension (4*N)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssyrfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO);
/**
*> \brief <b> SSYSV computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYSV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,<br>
*                         LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYSV computes the solution to a real system of linear equations<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          SSYTRF.<br>
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
*>          determined by SSYTRF.  If IPIV(k) > 0, then rows and columns<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >= 1, and for best performance<br>
*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for<br>
*>          SSYTRF.<br>
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
*> \ingroup realSYsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void ssysv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief <b> SSYSVX computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYSVX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysvx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysvx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysvx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,<br>
*                          LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK,<br>
*                          IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          FACT, UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS<br>
*       REAL               RCOND<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYSVX uses the diagonal pivoting factorization to compute the<br>
*> solution to a real system of linear equations A * X = B,<br>
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
*>          = 'F':  On entry, AF and IPIV contain the factored form of<br>
*>                  A.  AF and IPIV will not be modified.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>          If FACT = 'F', then AF is an input argument and on entry<br>
*>          contains the block diagonal matrix D and the multipliers used<br>
*>          to obtain the factor U or L from the factorization<br>
*>          A = U*D*U**T or A = L*D*L**T as computed by SSYTRF.<br>
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
*>          of D, as determined by SSYTRF.<br>
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
*>          of D, as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          The estimate of the reciprocal condition number of the matrix<br>
*>          A.  If RCOND is less than the machine precision (in<br>
*>          particular, if RCOND = 0), the matrix is singular to working<br>
*>          precision.  This condition is indicated by a return code of<br>
*>          INFO > 0.<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >= max(1,3*N), and for best<br>
*>          performance, when FACT = 'N', LWORK >= max(1,3*N,N*NB), where<br>
*>          NB is the optimal blocksize for SSYTRF.<br>
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
*> \ingroup realSYsolve<br>
*<br>
*  =====================================================================<br>
*/
	public void ssysvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
/**
*> \brief \b SSYSVXX<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYSVXX + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysvxx.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysvxx.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysvxx.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,<br>
*                           EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR,<br>
*                           N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP,<br>
*                           NPARAMS, PARAMS, WORK, IWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EQUED, FACT, UPLO<br>
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,<br>
*      $                   N_ERR_BNDS<br>
*       REAL               RCOND, RPVGRW<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * ), IWORK( * )<br>
*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),<br>
*      $                   X( LDX, * ), WORK( * )<br>
*       REAL               S( * ), PARAMS( * ), BERR( * ),<br>
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
*>    SSYSVXX uses the diagonal pivoting factorization to compute the<br>
*>    solution to a real system of linear equations A * X = B, where A<br>
*>    is an N-by-N symmetric matrix and X and B are N-by-NRHS matrices.<br>
*><br>
*>    If requested, both normwise and maximum componentwise error bounds<br>
*>    are returned. SSYSVXX will return a solution with a tiny<br>
*>    guaranteed error (O(eps) where eps is the working machine<br>
*>    precision) unless the matrix is very ill-conditioned, in which<br>
*>    case a warning is returned. Relevant condition numbers also are<br>
*>    calculated and returned.<br>
*><br>
*>    SSYSVXX accepts user-provided factorizations and equilibration<br>
*>    factors; see the definitions of the FACT and EQUED options.<br>
*>    Solving with refinement and using a factorization from a previous<br>
*>    SSYSVXX call will also produce a solution with either O(eps)<br>
*>    errors or warnings, but we cannot make that claim for general<br>
*>    user-provided factorizations and equilibration factors if they<br>
*>    differ from what SSYSVXX would itself produce.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          AF is REAL array, dimension (LDAF,N)<br>
*>     If FACT = 'F', then AF is an input argument and on entry<br>
*>     contains the block diagonal matrix D and the multipliers<br>
*>     used to obtain the factor U or L from the factorization A =<br>
*>     U*D*U**T or A = L*D*L**T as computed by SSYTRF.<br>
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
*>     structure of D, as determined by SSYTRF.  If IPIV(k) > 0,<br>
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
*>     structure of D, as determined by SSYTRF.<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          X is REAL array, dimension (LDX,NRHS)<br>
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
*>          WORK is REAL array, dimension (4*N)<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup realSYdriver<br>
*<br>
*  =====================================================================<br>
*/
	public void ssysvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL SLA_SYRPVGRW);
/**
*> \brief <b> SSYSV_ROOK computes the solution to system of linear equations A * X = B for SY matrices</b><br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYSV_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,<br>
*                         LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYSV_ROOK computes the solution to a real system of linear<br>
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
*> 1-by-1 and 2-by-2 diagonal blocks.<br>
*><br>
*> SSYTRF_ROOK is called to compute the factorization of a real<br>
*> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal<br>
*> pivoting method.<br>
*><br>
*> The factored form of A is then used to solve the system <br>
*> of equations A * X = B by calling SSYTRS_ROOK.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          SSYTRF_ROOK.<br>
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
*>          as determined by SSYTRF_ROOK.<br>
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
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The length of WORK.  LWORK >= 1, and for best performance<br>
*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for<br>
*>          SSYTRF_ROOK.<br>
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
*> \date April 2012<br>
*<br>
*> \ingroup realSYsolve<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   April 2012, Igor Kozachenko,<br>
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
	public void ssysv_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SSYSWAPR applies an elementary permutation on the rows and columns of a symmetric matrix.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYSWAPR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssyswapr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssyswapr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssyswapr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYSWAPR( UPLO, N, A, LDA, I1, I2)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER        UPLO<br>
*       INTEGER          I1, I2, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL             A( LDA, N )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYSWAPR applies an elementary permutation on the rows and the columns of<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the NB diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by SSYTRF.<br>
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
*> \ingroup realSYauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void ssyswapr_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER I1,INTEGER I2);
/**
*> \brief \b SSYTD2 reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTD2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytd2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytd2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytd2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), D( * ), E( * ), TAU( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal<br>
*> form T by an orthogonal similarity transformation: Q**T * A * Q = T.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          n-by-n upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading n-by-n lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal<br>
*>          of A are overwritten by the corresponding elements of the<br>
*>          tridiagonal matrix T, and the elements above the first<br>
*>          superdiagonal, with the array TAU, represent the orthogonal<br>
*>          matrix Q as a product of elementary reflectors; if UPLO<br>
*>          = 'L', the diagonal and first subdiagonal of A are over-<br>
*>          written by the corresponding elements of the tridiagonal<br>
*>          matrix T, and the elements below the first subdiagonal, with<br>
*>          the array TAU, represent the orthogonal matrix Q as a product<br>
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
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is REAL array, dimension (N-1)<br>
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
*> \ingroup realSYcomputational<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
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
	public void ssytd2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAU,INTEGER INFO);
/**
*> \brief \b SSYTF2 computes the factorization of a real symmetric indefinite matrix, using the diagonal pivoting method (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTF2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytf2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytf2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytf2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTF2( UPLO, N, A, LDA, IPIV, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTF2 computes the factorization of a real symmetric matrix A using<br>
*> the Bunch-Kaufman diagonal pivoting method:<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realSYcomputational<br>
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
*>    Replace l.204 and l.372<br>
*>         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN<br>
*>    by<br>
*>         IF( (MAX( ABSAKK, COLMAX ).EQ.ZERO) .OR. SISNAN(ABSAKK) ) THEN<br>
*><br>
*>  01-01-96 - Based on modifications by<br>
*>    J. Lewis, Boeing Computer Services Company<br>
*>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA<br>
*>  1-96 - Based on modifications by J. Lewis, Boeing Computer Services<br>
*>         Company<br>
*><br>
*> \endverbatim<br>
*<br>
*  =====================================================================<br>
*/
	public void ssytf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b SSYTF2_ROOK computes the factorization of a real symmetric indefinite matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method (unblocked algorithm).<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at<br>
*            http://www.netlib.org/lapack/explore-html/<br>
*<br>
*> \htmlonly<br>
*> Download SSYTF2_ROOK + dependencies<br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytf2_rook.f"><br>
*> [TGZ]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytf2_rook.f"><br>
*> [ZIP]</a><br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytf2_rook.f"><br>
*> [TXT]</a><br>
*> \endhtmlonly<br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTF2_ROOK( UPLO, N, A, LDA, IPIV, INFO )<br>
*<br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*<br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTF2_ROOK computes the factorization of a real symmetric matrix A<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*> \ingroup realSYcomputational<br>
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
	public void ssytf2_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
/**
*> \brief \b SSYTRD<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRD + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrd.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrd.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrd.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * ), D( * ), E( * ), TAU( * ),<br>
*      $                   WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRD reduces a real symmetric matrix A to real symmetric<br>
*> tridiagonal form T by an orthogonal similarity transformation:<br>
*> Q**T * A * Q = T.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading<br>
*>          N-by-N upper triangular part of A contains the upper<br>
*>          triangular part of the matrix A, and the strictly lower<br>
*>          triangular part of A is not referenced.  If UPLO = 'L', the<br>
*>          leading N-by-N lower triangular part of A contains the lower<br>
*>          triangular part of the matrix A, and the strictly upper<br>
*>          triangular part of A is not referenced.<br>
*>          On exit, if UPLO = 'U', the diagonal and first superdiagonal<br>
*>          of A are overwritten by the corresponding elements of the<br>
*>          tridiagonal matrix T, and the elements above the first<br>
*>          superdiagonal, with the array TAU, represent the orthogonal<br>
*>          matrix Q as a product of elementary reflectors; if UPLO<br>
*>          = 'L', the diagonal and first subdiagonal of A are over-<br>
*>          written by the corresponding elements of the tridiagonal<br>
*>          matrix T, and the elements below the first subdiagonal, with<br>
*>          the array TAU, represent the orthogonal matrix Q as a product<br>
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
*>          D is REAL array, dimension (N)<br>
*>          The diagonal elements of the tridiagonal matrix T:<br>
*>          D(i) = A(i,i).<br>
*> \endverbatim<br>
*><br>
*> \param[out] E<br>
*> \verbatim<br>
*>          E is REAL array, dimension (N-1)<br>
*>          The off-diagonal elements of the tridiagonal matrix T:<br>
*>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.<br>
*> \endverbatim<br>
*><br>
*> \param[out] TAU<br>
*> \verbatim<br>
*>          TAU is REAL array, dimension (N-1)<br>
*>          The scalar factors of the elementary reflectors (see Further<br>
*>          Details).<br>
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
*> \ingroup realSYcomputational<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
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
*>     H(i) = I - tau * v * v**T<br>
*><br>
*>  where tau is a real scalar, and v is a real vector with<br>
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
	public void ssytrd_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SSYTRF<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRF + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrf.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrf.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrf.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRF computes the factorization of a real symmetric matrix A using<br>
*> the Bunch-Kaufman diagonal pivoting method.  The form of the<br>
*> factorization is<br>
*><br>
*>    A = U*D*U**T  or  A = L*D*L**T<br>
*><br>
*> where U (or L) is a product of permutation and unit upper (lower)<br>
*> triangular matrices, and D is symmetric and block diagonal with <br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK))<br>
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
*> \ingroup realSYcomputational<br>
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
	public void ssytrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SSYTRF_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRF_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrf_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrf_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrf_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRF_ROOK computes the factorization of a real symmetric matrix A<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
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
*>          WORK is REAL array, dimension (MAX(1,LWORK)).<br>
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
*> \ingroup realSYcomputational<br>
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
	public void ssytrf_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SSYTRI<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRI + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRI computes the inverse of a real symmetric indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> SSYTRF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the block diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssytri_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSYTRI2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRI2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRI2 computes the inverse of a REAL symmetric indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> SSYTRF. SSYTRI2 sets the LEADING DIMENSION of the workspace<br>
*> before calling SSYTRI2X that actually computes the inverse.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the NB diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (N+NB+1)*(NB+3)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssytri2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b SSYTRI2X<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRI2X + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri2x.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri2x.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri2x.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRI2X( UPLO, N, A, LDA, IPIV, WORK, NB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N, NB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), WORK( N+NB+1,* )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRI2X computes the inverse of a real symmetric indefinite matrix<br>
*> A using the factorization A = U*D*U**T or A = L*D*L**T computed by<br>
*> SSYTRF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the NNB diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is REAL array, dimension (N+NNB+1,NNB+3)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssytri2x_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER NB,INTEGER INFO);
/**
*> \brief \b SSYTRI_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRI_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytri_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytri_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytri_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRI_ROOK computes the inverse of a real symmetric<br>
*> matrix A using the factorization A = U*D*U**T or A = L*D*L**T<br>
*> computed by SSYTRF_ROOK.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          On entry, the block diagonal matrix D and the multipliers<br>
*>          used to obtain the factor U or L as computed by SSYTRF_ROOK.<br>
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
*>          as determined by SSYTRF_ROOK.<br>
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
*> \date April 2012<br>
*<br>
*> \ingroup realSYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   April 2012, Igor Kozachenko,<br>
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
	public void ssytri_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSYTRS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRS solves a system of linear equations A*X = B with a real<br>
*> symmetric matrix A using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by SSYTRF.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssytrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
/**
*> \brief \b SSYTRS2<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRS2 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs2.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs2.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs2.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, <br>
*                           WORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), B( LDB, * ), WORK( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRS2 solves a system of linear equations A*X = B with a real<br>
*> symmetric matrix A using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by SSYTRF and converted by SSYCONV.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSYTRF.<br>
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
*>          as determined by SSYTRF.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*>          WORK is REAL array, dimension (N)<br>
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
*> \ingroup realSYcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void ssytrs2_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,float[] WORK,INTEGER INFO);
/**
*> \brief \b SSYTRS_ROOK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SSYTRS_ROOK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssytrs_rook.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssytrs_rook.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssytrs_rook.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       INTEGER            INFO, LDA, LDB, N, NRHS<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       INTEGER            IPIV( * )<br>
*       REAL               A( LDA, * ), B( LDB, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SSYTRS_ROOK solves a system of linear equations A*X = B with<br>
*> a real symmetric matrix A using the factorization A = U*D*U**T or<br>
*> A = L*D*L**T computed by SSYTRF_ROOK.<br>
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
*>          A is REAL array, dimension (LDA,N)<br>
*>          The block diagonal matrix D and the multipliers used to<br>
*>          obtain the factor U or L as computed by SSYTRF_ROOK.<br>
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
*>          as determined by SSYTRF_ROOK.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] B<br>
*> \verbatim<br>
*>          B is REAL array, dimension (LDB,NRHS)<br>
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
*> \date April 2012<br>
*<br>
*> \ingroup realSYcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> \verbatim<br>
*><br>
*>   April 2012, Igor Kozachenko,<br>
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
	public void ssytrs_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);

}