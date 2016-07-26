package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackLS extends Library
{

	public static LapackLS instance = (LapackLS) Native.loadLibrary("liblapack",LapackLS.class);

/**
*> \brief \b LSAMEN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download LSAMEN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/lsamen.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/lsamen.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/lsamen.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       LOGICAL          FUNCTION LSAMEN( N, CA, CB )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*( * )    CA, CB<br>
*       INTEGER            N<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> LSAMEN  tests if the first N letters of CA are the same as the<br>
*> first N letters of CB, regardless of case.<br>
*> LSAMEN returns .TRUE. if CA and CB are equivalent except for case<br>
*> and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )<br>
*> or LEN( CB ) is less than N.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of characters in CA and CB to be compared.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CA<br>
*> \verbatim<br>
*>          CA is CHARACTER*(*)<br>
*> \endverbatim<br>
*><br>
*> \param[in] CB<br>
*> \verbatim<br>
*>          CB is CHARACTER*(*)<br>
*>          CA and CB specify two character strings of length at least N.<br>
*>          Only the first N characters of each string will be accessed.<br>
*> \endverbatim<br>
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
	public boolean lsamen_(INTEGER N,char[] CA,char[] CB);

}