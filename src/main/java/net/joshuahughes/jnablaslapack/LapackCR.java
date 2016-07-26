package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCR extends Library
{

	public static LapackCR instance = (LapackCR) Native.loadLibrary("liblapack",LapackCR.class);

/**
*> \brief \b CROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download CROT + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/crot.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/crot.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/crot.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            INCX, INCY, N<br>
*       REAL               C<br>
*       COMPLEX            S<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            CX( * ), CY( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> CROT   applies a plane rotation, where the cos (C) is real and the<br>
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
*>          CX is COMPLEX array, dimension (N)<br>
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
*>          CY is COMPLEX array, dimension (N)<br>
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
*>          C is REAL<br>
*> \endverbatim<br>
*><br>
*> \param[in] S<br>
*> \verbatim<br>
*>          S is COMPLEX<br>
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
*> \ingroup complexOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void crot_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY,REAL C,float[] S);

}