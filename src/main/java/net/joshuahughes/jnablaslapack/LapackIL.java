package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIL extends Library
{

	public static LapackIL instance = (LapackIL) Native.loadLibrary("liblapack",LapackIL.class);

/**
*> \brief \b ILACLC scans a matrix for its last non-zero column.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILACLC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaclc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaclc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaclc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILACLC( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILACLC scans A for its last non-zero column.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int ilaclc_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
/**
*> \brief \b ILACLR scans a matrix for its last non-zero row.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILACLR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaclr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaclr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaclr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILACLR( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX            A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILACLR scans A for its last non-zero row.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int ilaclr_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
/**
*> \brief \b ILADIAG<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILADIAG + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladiag.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladiag.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladiag.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILADIAG( DIAG )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          DIAG<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine translated from a character string specifying if a<br>
*> matrix has unit diagonal or not to the relevant BLAST-specified<br>
*> integer constant.<br>
*><br>
*> ILADIAG returns an INTEGER.  If ILADIAG < 0, then the input is not a<br>
*> character indicating a unit or non-unit diagonal.  Otherwise ILADIAG<br>
*> returns the constant value corresponding to DIAG.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public int iladiag_(CHARACTER DIAG);
/**
*> \brief \b ILADLC scans a matrix for its last non-zero column.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILADLC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILADLC( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILADLC scans A for its last non-zero column.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int iladlc_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);
/**
*> \brief \b ILADLR scans a matrix for its last non-zero row.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILADLR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladlr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladlr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILADLR( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILADLR scans A for its last non-zero row.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is DOUBLE PRECISION array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int iladlr_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);
/**
*> \brief \b ILAENV<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILAENV + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER*( * )    NAME, OPTS<br>
*       INTEGER            ISPEC, N1, N2, N3, N4<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILAENV is called from the LAPACK routines to choose problem-dependent<br>
*> parameters for the local environment.  See ISPEC for a description of<br>
*> the parameters.<br>
*><br>
*> ILAENV returns an INTEGER<br>
*> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC<br>
*> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.<br>
*><br>
*> This version provides a set of parameters which should give good,<br>
*> but not optimal, performance on many of the currently available<br>
*> computers.  Users are encouraged to modify this subroutine to set<br>
*> the tuning parameters for their particular machine using the option<br>
*> and problem size information in the arguments.<br>
*><br>
*> This routine will not function correctly if it is converted to all<br>
*> lower case.  Converting it to all upper case is allowed.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] ISPEC<br>
*> \verbatim<br>
*>          ISPEC is INTEGER<br>
*>          Specifies the parameter to be returned as the value of<br>
*>          ILAENV.<br>
*>          = 1: the optimal blocksize; if this value is 1, an unblocked<br>
*>               algorithm will give the best performance.<br>
*>          = 2: the minimum block size for which the block routine<br>
*>               should be used; if the usable block size is less than<br>
*>               this value, an unblocked routine should be used.<br>
*>          = 3: the crossover point (in a block routine, for N less<br>
*>               than this value, an unblocked routine should be used)<br>
*>          = 4: the number of shifts, used in the nonsymmetric<br>
*>               eigenvalue routines (DEPRECATED)<br>
*>          = 5: the minimum column dimension for blocking to be used;<br>
*>               rectangular blocks must have dimension at least k by m,<br>
*>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)<br>
*>          = 6: the crossover point for the SVD (when reducing an m by n<br>
*>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds<br>
*>               this value, a QR factorization is used first to reduce<br>
*>               the matrix to a triangular form.)<br>
*>          = 7: the number of processors<br>
*>          = 8: the crossover point for the multishift QR method<br>
*>               for nonsymmetric eigenvalue problems (DEPRECATED)<br>
*>          = 9: maximum size of the subproblems at the bottom of the<br>
*>               computation tree in the divide-and-conquer algorithm<br>
*>               (used by xGELSD and xGESDD)<br>
*>          =10: ieee NaN arithmetic can be trusted not to trap<br>
*>          =11: infinity arithmetic can be trusted not to trap<br>
*>          12 <= ISPEC <= 16:<br>
*>               xHSEQR or related subroutines,<br>
*>               see IPARMQ for detailed explanation<br>
*> \endverbatim<br>
*><br>
*> \param[in] NAME<br>
*> \verbatim<br>
*>          NAME is CHARACTER*(*)<br>
*>          The name of the calling subroutine, in either upper case or<br>
*>          lower case.<br>
*> \endverbatim<br>
*><br>
*> \param[in] OPTS<br>
*> \verbatim<br>
*>          OPTS is CHARACTER*(*)<br>
*>          The character options to the subroutine NAME, concatenated<br>
*>          into a single character string.  For example, UPLO = 'U',<br>
*>          TRANS = 'T', and DIAG = 'N' for a triangular routine would<br>
*>          be specified as OPTS = 'UTN'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N1<br>
*> \verbatim<br>
*>          N1 is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] N2<br>
*> \verbatim<br>
*>          N2 is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] N3<br>
*> \verbatim<br>
*>          N3 is INTEGER<br>
*> \endverbatim<br>
*><br>
*> \param[in] N4<br>
*> \verbatim<br>
*>          N4 is INTEGER<br>
*>          Problem dimensions for the subroutine NAME; these may not all<br>
*>          be required.<br>
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
*> \date June 2016<br>
*<br>
*> \ingroup auxOTHERauxiliary<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  The following conventions have been used when calling ILAENV from the<br>
*>  LAPACK routines:<br>
*>  1)  OPTS is a concatenation of all of the character options to<br>
*>      subroutine NAME, in the same order that they appear in the<br>
*>      argument list for NAME, even if they are not used in determining<br>
*>      the value of the parameter specified by ISPEC.<br>
*>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order<br>
*>      that they appear in the argument list for NAME.  N1 is used<br>
*>      first, N2 second, and so on, and unused problem dimensions are<br>
*>      passed a value of -1.<br>
*>  3)  The parameter value returned by ILAENV is checked for validity in<br>
*>      the calling subroutine.  For example, ILAENV is used to retrieve<br>
*>      the optimal blocksize for STRTRI as follows:<br>
*><br>
*>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )<br>
*>      IF( NB.LE.1 ) NB = MAX( 1, N )<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public int ilaenv_(INTEGER ISPEC,byte[] NAME,byte[] OPTS,INTEGER N1,INTEGER N2,INTEGER N3,INTEGER N4);
/**
*> \brief \b ILAPREC<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILAPREC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaprec.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaprec.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaprec.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILAPREC( PREC )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          PREC<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine translated from a character string specifying an<br>
*> intermediate precision to the relevant BLAST-specified integer<br>
*> constant.<br>
*><br>
*> ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a<br>
*> character indicating a supported intermediate precision.  Otherwise<br>
*> ILAPREC returns the constant value corresponding to PREC.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public int ilaprec_(CHARACTER PREC);
/**
*> \brief \b ILASLC scans a matrix for its last non-zero column.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILASLC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILASLC( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILASLC scans A for its last non-zero column.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int ilaslc_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
/**
*> \brief \b ILASLR scans a matrix for its last non-zero row.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILASLR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILASLR( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       REAL               A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILASLR scans A for its last non-zero row.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is REAL array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int ilaslr_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
/**
*> \brief \b ILATRANS<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILATRANS + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilatrans.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilatrans.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilatrans.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILATRANS( TRANS )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          TRANS<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine translates from a character string specifying a<br>
*> transposition operation to the relevant BLAST-specified integer<br>
*> constant.<br>
*><br>
*> ILATRANS returns an INTEGER.  If ILATRANS < 0, then the input is not<br>
*> a character indicating a transposition operator.  Otherwise ILATRANS<br>
*> returns the constant value corresponding to TRANS.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public int ilatrans_(CHARACTER TRANS);
/**
*> \brief \b ILAUPLO<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILAUPLO + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilauplo.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilauplo.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilauplo.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILAUPLO( UPLO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          UPLO<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> This subroutine translated from a character string specifying a<br>
*> upper- or lower-triangular matrix to the relevant BLAST-specified<br>
*> integer constant.<br>
*><br>
*> ILAUPLO returns an INTEGER.  If ILAUPLO < 0, then the input is not<br>
*> a character indicating an upper- or lower-triangular matrix.<br>
*> Otherwise ILAUPLO returns the constant value corresponding to UPLO.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
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
*> \ingroup auxOTHERcomputational<br>
*<br>
*  =====================================================================<br>
*/
	public int ilauplo_(CHARACTER UPLO);
/**
*> \brief \b ILAVER returns the LAPACK version.<br>
**<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*     SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )<br>
*<br>
*     INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>  This subroutine returns the LAPACK version.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*>  \param[out] VERS_MAJOR<br>
*>  \verbatim<br>
*>      return the lapack major version<br>
*>  \endverbatim<br>
*><br>
*>  \param[out] VERS_MINOR<br>
*>  \verbatim<br>
*>      return the lapack minor version from the major version<br>
*>  \endverbatim<br>
*><br>
*>  \param[out] VERS_PATCH<br>
*>  \verbatim<br>
*>      return the lapack patch version from the minor version<br>
*>  \endverbatim<br>
*<br>
*  Authors:<br>
*  ========<br>
*<br>
*> \author Univ. of Tennessee <br>
*> \author Univ. of California Berkeley <br>
*> \author Univ. of Colorado Denver <br>
*> \author NAG Ltd. <br>
*<br>
*> \date June 2016<br>
*<br>
*> \ingroup auxOTHERauxiliary<br>
*<br>
*  =====================================================================<br>
*/
	public void ilaver_(INTEGER VERS_MAJOR,INTEGER VERS_MINOR,INTEGER VERS_PATCH);
/**
*> \brief \b ILAZLC scans a matrix for its last non-zero column.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILAZLC + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlc.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlc.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlc.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILAZLC( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILAZLC scans A for its last non-zero column.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int ilazlc_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);
/**
*> \brief \b ILAZLR scans a matrix for its last non-zero row.<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download ILAZLR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       INTEGER FUNCTION ILAZLR( M, N, A, LDA )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            M, N, LDA<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       COMPLEX*16         A( LDA, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> ILAZLR scans A for its last non-zero row.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of rows of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The number of columns of the matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] A<br>
*> \verbatim<br>
*>          A is COMPLEX*16 array, dimension (LDA,N)<br>
*>          The m by n matrix A.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDA<br>
*> \verbatim<br>
*>          LDA is INTEGER<br>
*>          The leading dimension of the array A. LDA >= max(1,M).<br>
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
	public int ilazlr_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);

}