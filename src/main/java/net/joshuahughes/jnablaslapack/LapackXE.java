package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackXE extends Library
{

	public static LapackXE instance = (LapackXE) Native.loadLibrary("liblapack",LapackXE.class);

/**
*> \brief \b XERBLA<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download XERBLA + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE XERBLA( SRNAME, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*(*)      SRNAME<br>
*       INTEGER            INFO<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> XERBLA  is an error handler for the LAPACK routines.<br>
*> It is called by an LAPACK routine if an input parameter has an<br>
*> invalid value.  A message is printed and execution stops.<br>
*><br>
*> Installers may consider modifying the STOP statement in order to<br>
*> call system-specific exception-handling facilities.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SRNAME<br>
*> \verbatim<br>
*>          SRNAME is CHARACTER*(*)<br>
*>          The name of the routine which called XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          The position of the invalid parameter in the parameter list<br>
*>          of the calling routine.<br>
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
	public void xerbla_(char[] SRNAME,INTEGER INFO);
/**
*> \brief \b XERBLA_ARRAY<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download XERBLA_ARRAY + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla_array.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla_array.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla_array.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE XERBLA_ARRAY( SRNAME_ARRAY, SRNAME_LEN, INFO)<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER SRNAME_LEN, INFO<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       CHARACTER(1) SRNAME_ARRAY(SRNAME_LEN)<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> XERBLA_ARRAY assists other languages in calling XERBLA, the LAPACK<br>
*> and BLAS error handler.  Rather than taking a Fortran string argument<br>
*> as the function's name, XERBLA_ARRAY takes an array of single<br>
*> characters along with the array's length.  XERBLA_ARRAY then copies<br>
*> up to 32 characters of that array into a Fortran string and passes<br>
*> that to XERBLA.  If called with a non-positive SRNAME_LEN,<br>
*> XERBLA_ARRAY will call XERBLA with a string of all blank characters.<br>
*><br>
*> Say some macro or other device makes XERBLA_ARRAY available to C99<br>
*> by a name lapack_xerbla and with a common Fortran calling convention.<br>
*> Then a C99 program could invoke XERBLA via:<br>
*>    {<br>
*>      int flen = strlen(__func__);<br>
*>      lapack_xerbla(__func__, &flen, &info);<br>
*>    }<br>
*><br>
*> Providing XERBLA_ARRAY is not necessary for intercepting LAPACK<br>
*> errors.  XERBLA_ARRAY calls XERBLA.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SRNAME_ARRAY<br>
*> \verbatim<br>
*>          SRNAME_ARRAY is CHARACTER(1) array, dimension (SRNAME_LEN)<br>
*>          The name of the routine which called XERBLA_ARRAY.<br>
*> \endverbatim<br>
*><br>
*> \param[in] SRNAME_LEN<br>
*> \verbatim<br>
*>          SRNAME_LEN is INTEGER<br>
*>          The length of the name in SRNAME_ARRAY.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          The position of the invalid parameter in the parameter list<br>
*>          of the calling routine.<br>
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
	public void xerbla_array_(char[] SRNAME_ARRAY,INTEGER SRNAME_LEN,INTEGER INFO);

}