package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSR extends Library
{

	public static LapackSR instance = (LapackSR) Native.loadLibrary("liblapack",LapackSR.class);

/**
*> \brief \b SRSCL multiplies a vector by the reciprocal of a real scalar.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download SRSCL + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/srscl.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/srscl.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/srscl.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE SRSCL( N, SA, SX, INCX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       REAL               SA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               SX( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> SRSCL multiplies an n-element real vector x by the real scalar 1/a.<br>
*> This is done without overflow or underflow as long as<br>
*> the final result x/a does not overflow or underflow.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of components of the vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SA<br>
*> \verbatim<br>
*>          SA is REAL<br>
*>          The scalar a which is used to divide each component of x.<br>
*>          SA must be >= 0, or the subroutine will divide by zero.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SX<br>
*> \verbatim<br>
*>          SX is REAL array, dimension<br>
*>                         (1+(N-1)*abs(INCX))<br>
*>          The n-element vector x.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between successive values of the vector SX.<br>
*>          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n<br>
*> \endverbatim<br>
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
*  =====================================================================<br>
*/
	public void srscl_(INTEGER N,REAL SA,float[] SX,INTEGER INCX);

}