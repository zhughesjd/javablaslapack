package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIP extends Library
{

	public static LapackIP instance = (LapackIP) Native.loadLibrary("liblapack",LapackIP.class);

/**
*> \brief \b IPARMQ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download IPARMQ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, ISPEC, LWORK, N<br>
*       CHARACTER          NAME*( * ), OPTS*( * )<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>      This program sets problem and machine dependent parameters<br>
*>      useful for xHSEQR and related subroutines for eigenvalue<br>
*>      problems. It is called whenever<br>
*>      IPARMQ is called with 12 <= ISPEC <= 16<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ISPEC<br>
*> \verbatim<br>
*>          ISPEC is integer scalar<br>
*>              ISPEC specifies which tunable parameter IPARMQ should<br>
*>              return.<br>
*><br>
*>              ISPEC=12: (INMIN)  Matrices of order nmin or less<br>
*>                        are sent directly to xLAHQR, the implicit<br>
*>                        double shift QR algorithm.  NMIN must be<br>
*>                        at least 11.<br>
*><br>
*>              ISPEC=13: (INWIN)  Size of the deflation window.<br>
*>                        This is best set greater than or equal to<br>
*>                        the number of simultaneous shifts NS.<br>
*>                        Larger matrices benefit from larger deflation<br>
*>                        windows.<br>
*><br>
*>              ISPEC=14: (INIBL) Determines when to stop nibbling and<br>
*>                        invest in an (expensive) multi-shift QR sweep.<br>
*>                        If the aggressive early deflation subroutine<br>
*>                        finds LD converged eigenvalues from an order<br>
*>                        NW deflation window and LD.GT.(NW*NIBBLE)/100,<br>
*>                        then the next QR sweep is skipped and early<br>
*>                        deflation is applied immediately to the<br>
*>                        remaining active diagonal block.  Setting<br>
*>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a<br>
*>                        multi-shift QR sweep whenever early deflation<br>
*>                        finds a converged eigenvalue.  Setting<br>
*>                        IPARMQ(ISPEC=14) greater than or equal to 100<br>
*>                        prevents TTQRE from skipping a multi-shift<br>
*>                        QR sweep.<br>
*><br>
*>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in<br>
*>                        a multi-shift QR iteration.<br>
*><br>
*>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the<br>
*>                        following meanings.<br>
*>                        0:  During the multi-shift QR/QZ sweep,<br>
*>                            blocked eigenvalue reordering, blocked<br>
*>                            Hessenberg-triangular reduction,<br>
*>                            reflections and/or rotations are not<br>
*>                            accumulated when updating the<br>
*>                            far-from-diagonal matrix entries.<br>
*>                        1:  During the multi-shift QR/QZ sweep,<br>
*>                            blocked eigenvalue reordering, blocked<br>
*>                            Hessenberg-triangular reduction,<br>
*>                            reflections and/or rotations are<br>
*>                            accumulated, and matrix-matrix<br>
*>                            multiplication is used to update the<br>
*>                            far-from-diagonal matrix entries.<br>
*>                        2:  During the multi-shift QR/QZ sweep,<br>
*>                            blocked eigenvalue reordering, blocked<br>
*>                            Hessenberg-triangular reduction,<br>
*>                            reflections and/or rotations are<br>
*>                            accumulated, and 2-by-2 block structure<br>
*>                            is exploited during matrix-matrix<br>
*>                            multiplies.<br>
*>                        (If xTRMM is slower than xGEMM, then<br>
*>                        IPARMQ(ISPEC=16)=1 may be more efficient than<br>
*>                        IPARMQ(ISPEC=16)=2 despite the greater level of<br>
*>                        arithmetic work implied by the latter choice.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] NAME<br>
*> \verbatim<br>
*>          NAME is character string<br>
*>               Name of the calling subroutine<br>
*> \endverbatim<br>
*><br>
*> \param[in] OPTS<br>
*> \verbatim<br>
*>          OPTS is character string<br>
*>               This is a concatenation of the string arguments to<br>
*>               TTQRE.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is integer scalar<br>
*>               N is the order of the Hessenberg matrix H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] ILO<br>
*> \verbatim<br>
*>          ILO is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] IHI<br>
*> \verbatim<br>
*>          IHI is INTEGER<br>
*>               It is assumed that H is already upper triangular<br>
*>               in rows and columns 1:ILO-1 and IHI+1:N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is integer scalar<br>
*>               The amount of workspace available.<br>
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
*> \date November 2015<br>
*<br>
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>       Little is known about how best to choose these parameters.<br>
*>       It is possible to use different values of the parameters<br>
*>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.<br>
*><br>
*>       It is probably best to choose different parameters for<br>
*>       different matrices and different parameters at different<br>
*>       times during the iteration, but this has not been<br>
*>       implemented --- yet.<br>
*><br>
*><br>
*>       The best choices of most of the parameters depend<br>
*>       in an ill-understood way on the relative execution<br>
*>       rate of xLAQR3 and xLAQR5 and on the nature of each<br>
*>       particular eigenvalue problem.  Experiment may be the<br>
*>       only practical way to determine which choices are most<br>
*>       effective.<br>
*><br>
*>       Following is a list of default values supplied by IPARMQ.<br>
*>       These defaults may be adjusted in order to attain better<br>
*>       performance in any particular computational environment.<br>
*><br>
*>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.<br>
*>                        Default: 75. (Must be at least 11.)<br>
*><br>
*>       IPARMQ(ISPEC=13) Recommended deflation window size.<br>
*>                        This depends on ILO, IHI and NS, the<br>
*>                        number of simultaneous shifts returned<br>
*>                        by IPARMQ(ISPEC=15).  The default for<br>
*>                        (IHI-ILO+1).LE.500 is NS.  The default<br>
*>                        for (IHI-ILO+1).GT.500 is 3*NS/2.<br>
*><br>
*>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.<br>
*><br>
*>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.<br>
*>                        a multi-shift QR iteration.<br>
*><br>
*>                        If IHI-ILO+1 is ...<br>
*><br>
*>                        greater than      ...but less    ... the<br>
*>                        or equal to ...      than        default is<br>
*><br>
*>                                0               30       NS =   2+<br>
*>                               30               60       NS =   4+<br>
*>                               60              150       NS =  10<br>
*>                              150              590       NS =  **<br>
*>                              590             3000       NS =  64<br>
*>                             3000             6000       NS = 128<br>
*>                             6000             infinity   NS = 256<br>
*><br>
*>                    (+)  By default matrices of this order are<br>
*>                         passed to the implicit double shift routine<br>
*>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These<br>
*>                         values of NS are used only in case of a rare<br>
*>                         xLAHQR failure.<br>
*><br>
*>                    (**) The asterisks (**) indicate an ad-hoc<br>
*>                         function increasing from 10 to 64.<br>
*><br>
*>       IPARMQ(ISPEC=16) Select structured matrix multiply.<br>
*>                        (See ISPEC=16 above for details.)<br>
*>                        Default: 3.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public int iparmq_(INTEGER ISPEC,CHARACTER NAME,CHARACTER OPTS,INTEGER N,INTEGER ILO,INTEGER IHI,INTEGER LWORK);

}