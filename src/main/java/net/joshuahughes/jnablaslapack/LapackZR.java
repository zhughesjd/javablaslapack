package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZR extends Library
{

	public static LapackZR instance = (LapackZR) Native.loadLibrary("liblapack",LapackZR.class);

/**
*> \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ZROT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, INCY, N<br>
*       DOUBLE PRECISION   C<br>
*       COMPLEX*16         S<br>
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
*> ZROT   applies a plane rotation, where the cos (C) is real and the<br>
*> sin (S) is complex, and the vectors CX and CY are complex.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of elements in the vectors CX and CY.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CX<br>
*> \verbatim<br>
*>          CX is COMPLEX*16 array, dimension (N)<br>
*>          On input, the vector X.<br>
*>          On output, CX is overwritten with C*X + S*Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCX<br>
*> \verbatim<br>
*>          INCX is INTEGER<br>
*>          The increment between successive values of CY.  INCX <> 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] CY<br>
*> \verbatim<br>
*>          CY is COMPLEX*16 array, dimension (N)<br>
*>          On input, the vector Y.<br>
*>          On output, CY is overwritten with -CONJG(S)*X + C*Y.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INCY<br>
*> \verbatim<br>
*>          INCY is INTEGER<br>
*>          The increment between successive values of CY.  INCX <> 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] C<br>
*> \verbatim<br>
*>          C is DOUBLE PRECISION<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is COMPLEX*16<br>
*>          C and S define a rotation<br>
*>             [  C          S  ]<br>
*>             [ -conjg(S)   C  ]<br>
*>          where C*C + S*CONJG(S) = 1.0.<br>
*> \endverbatim<br>
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
*  =====================================================================<br>
*/
	public void zrot_(INTEGER N,double[] CX,INTEGER INCX,double[] CY,INTEGER INCY,DOUBLE C,double[] S);

}