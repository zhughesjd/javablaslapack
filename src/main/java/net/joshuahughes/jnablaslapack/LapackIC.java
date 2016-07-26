package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIC extends Library
{

	public static LapackIC instance = (LapackIC) Native.loadLibrary("liblapack",LapackIC.class);

/**
*> \brief \b ICMAX1 finds the index of the first vector element of maximum absolute value.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ICMAX1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/icmax1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/icmax1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/icmax1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER          FUNCTION ICMAX1( N, CX, INCX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            CX( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ICMAX1 finds the index of the first vector element of maximum absolute value.<br>
*><br>
*> Based on ICAMAX from Level 1 BLAS.<br>
*> The change is to use the 'genuine' absolute value.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of elements in the vector CX.<br>
*> \endverbatim<br>
*><br>
*> \param[in] CX<br>
*> \verbatim<br>
*>          CX is COMPLEX array, dimension (N)<br>
*>          The vector CX. The ICMAX1 function returns the index of its first<br>
*>          element of maximum absolute value.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The spacing between successive values of CX.  INCX >= 1.<br>
*> \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date February 2014<br>
*<br>
*> \ingroup complexOTHERauxiliary<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Nick Higham for use with CLACON.<br>
*<br>
*  =====================================================================<br>
*/
	public int icmax1_(INTEGER N,float[] CX,INTEGER INCX);

}