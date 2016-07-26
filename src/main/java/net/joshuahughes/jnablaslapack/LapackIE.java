package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIE extends Library
{

	public static LapackIE instance = (LapackIE) Native.loadLibrary("liblapack",LapackIE.class);

/**
*> \brief \b IEEECK<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download IEEECK + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            ISPEC<br>
*       REAL               ONE, ZERO<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> IEEECK is called from the ILAENV to verify that Infinity and<br>
*> possibly NaN arithmetic is safe (i.e. will not trap).<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ISPEC<br>
*> \verbatim<br>
*>          ISPEC is INTEGER<br>
*>          Specifies whether to test just for inifinity arithmetic<br>
*>          or whether to test for infinity and NaN arithmetic.<br>
*>          = 0: Verify infinity arithmetic only.<br>
*>          = 1: Verify infinity and NaN arithmetic.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ZERO<br>
*> \verbatim<br>
*>          ZERO is REAL<br>
*>          Must contain the value 0.0<br>
*>          This is passed to prevent the compiler from optimizing<br>
*>          away this code.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ONE<br>
*> \verbatim<br>
*>          ONE is REAL<br>
*>          Must contain the value 1.0<br>
*>          This is passed to prevent the compiler from optimizing<br>
*>          away this code.<br>
*><br>
*>  RETURN VALUE:  INTEGER<br>
*>          = 0:  Arithmetic failed to produce the correct answers<br>
*>          = 1:  Arithmetic produced the correct answers<br>
*> \endverbatim<br>
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
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public int ieeeck_(INTEGER ISPEC,REAL ZERO,REAL ONE);

}