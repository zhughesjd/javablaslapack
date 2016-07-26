package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDI extends Library
{

	public static LapackDI instance = (LapackDI) Native.loadLibrary("liblapack",LapackDI.class);

/**
*> \brief \b DISNAN tests input for NaN.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DISNAN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       LOGICAL FUNCTION DISNAN( DIN )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       DOUBLE PRECISION   DIN<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.<br>
*> otherwise.  To be replaced by the Fortran 2003 intrinsic in the<br>
*> future.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] DIN<br>
*> \verbatim<br>
*>          DIN is DOUBLE PRECISION<br>
*>          Input to test for NaN.<br>
*> \endverbatim<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public boolean disnan_(DOUBLE DIN);

}