package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSD extends Library
{

	public static LapackSD instance = (LapackSD) Native.loadLibrary("liblapack",LapackSD.class);

/**
*> \brief \b SDISNA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SDISNA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sdisna.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sdisna.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sdisna.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SDISNA( JOB, M, N, D, SEP, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          JOB<br>
*       INTEGER            INFO, M, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               D( * ), SEP( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SDISNA computes the reciprocal condition numbers for the eigenvectors<br>
*> of a real symmetric or complex Hermitian matrix or for the left or<br>
*> right singular vectors of a general m-by-n matrix. The reciprocal<br>
*> condition number is the 'gap' between the corresponding eigenvalue or<br>
*> singular value and the nearest other one.<br>
*><br>
*> The bound on the error, measured by angle in radians, in the I-th<br>
*> computed vector is given by<br>
*><br>
*>        SLAMCH( 'E' ) * ( ANORM / SEP( I ) )<br>
*><br>
*> where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed<br>
*> to be smaller than SLAMCH( 'E' )*ANORM in order to limit the size of<br>
*> the error bound.<br>
*><br>
*> SDISNA may also be used to compute error bounds for eigenvectors of<br>
*> the generalized symmetric definite eigenproblem.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          Specifies for which problem the reciprocal condition numbers<br>
*>          should be computed:<br>
*>          = 'E':  the eigenvectors of a symmetric/Hermitian matrix;<br>
*>          = 'L':  the left singular vectors of a general matrix;<br>
*>          = 'R':  the right singular vectors of a general matrix.<br>
*> \endverbatim<br>
*><br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix. M >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          If JOB = 'L' or 'R', the number of columns of the matrix,<br>
*>          in which case N >= 0. Ignored if JOB = 'E'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] D<br>
*> \verbatim<br>
*>          D is REAL array, dimension (M) if JOB = 'E'<br>
*>                              dimension (min(M,N)) if JOB = 'L' or 'R'<br>
*>          The eigenvalues (if JOB = 'E') or singular values (if JOB =<br>
*>          'L' or 'R') of the matrix, in either increasing or decreasing<br>
*>          order. If singular values, they must be non-negative.<br>
*> \endverbatim<br>
*><br>
*> \param[out] SEP<br>
*> \verbatim<br>
*>          SEP is REAL array, dimension (M) if JOB = 'E'<br>
*>                               dimension (min(M,N)) if JOB = 'L' or 'R'<br>
*>          The reciprocal condition numbers of the vectors.<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public void sdisna_(CHARACTER JOB,INTEGER M,INTEGER N,float[] D,float[] SEP,INTEGER INFO);

}