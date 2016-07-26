package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDZ extends Library
{

	public static LapackDZ instance = (LapackDZ) Native.loadLibrary("liblapack",LapackDZ.class);

/**
*> \brief \b DZSUM1 forms the 1-norm of the complex vector using the true absolute value.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DZSUM1 + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dzsum1.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dzsum1.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dzsum1.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       DOUBLE PRECISION FUNCTION DZSUM1( N, CX, INCX )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         CX( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DZSUM1 takes the sum of the absolute values of a complex<br>
*> vector and returns a double precision result.<br>
*><br>
*> Based on DZASUM from the Level 1 BLAS.<br>
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
*>          CX is COMPLEX*16 array, dimension (N)<br>
*>          The vector whose elements will be summed.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The spacing between successive values of CX.  INCX > 0.<br>
*> \endverbatim<br>
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
*> \par Contributors:<br>
*  ==================<br>
*><br>
*> Nick Higham for use with ZLACON.<br>
*<br>
*  =====================================================================<br>
*/
	public double dzsum1_(INTEGER N,double[] CX,INTEGER INCX);

}