package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDH extends Library
{

	public static LapackDH instance = (LapackDH) Native.loadLibrary("liblapack",LapackDH.class);

/**
*> \brief \b DHGEQZ<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DHGEQZ + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhgeqz.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhgeqz.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhgeqz.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,<br>
*                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,<br>
*                          LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          COMPQ, COMPZ, JOB<br>
*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ),<br>
*      $                   H( LDH, * ), Q( LDQ, * ), T( LDT, * ),<br>
*      $                   WORK( * ), Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DHGEQZ computes the eigenvalues of a real matrix pair (H,T),<br>
*> where H is an upper Hessenberg matrix and T is upper triangular,<br>
*> using the double-shift QZ method.<br>
*> Matrix pairs of this type are produced by the reduction to<br>
*> generalized upper Hessenberg form of a real matrix pair (A,B):<br>
*><br>
*>    A = Q1*H*Z1**T,  B = Q1*T*Z1**T,<br>
*><br>
*> as computed by DGGHRD.<br>
*><br>
*> If JOB='S', then the Hessenberg-triangular pair (H,T) is<br>
*> also reduced to generalized Schur form,<br>
*> <br>
*>    H = Q*S*Z**T,  T = Q*P*Z**T,<br>
*> <br>
*> where Q and Z are orthogonal matrices, P is an upper triangular<br>
*> matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2<br>
*> diagonal blocks.<br>
*><br>
*> The 1-by-1 blocks correspond to real eigenvalues of the matrix pair<br>
*> (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of<br>
*> eigenvalues.<br>
*><br>
*> Additionally, the 2-by-2 upper triangular diagonal blocks of P<br>
*> corresponding to 2-by-2 blocks of S are reduced to positive diagonal<br>
*> form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,<br>
*> P(j,j) > 0, and P(j+1,j+1) > 0.<br>
*><br>
*> Optionally, the orthogonal matrix Q from the generalized Schur<br>
*> factorization may be postmultiplied into an input matrix Q1, and the<br>
*> orthogonal matrix Z may be postmultiplied into an input matrix Z1.<br>
*> If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced<br>
*> the matrix pair (A,B) to generalized upper Hessenberg form, then the<br>
*> output matrices Q1*Q and Z1*Z are the orthogonal factors from the<br>
*> generalized Schur factorization of (A,B):<br>
*><br>
*>    A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.<br>
*> <br>
*> To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,<br>
*> of (A,B)) are computed as a pair of values (alpha,beta), where alpha is<br>
*> complex and beta real.<br>
*> If beta is nonzero, lambda = alpha / beta is an eigenvalue of the<br>
*> generalized nonsymmetric eigenvalue problem (GNEP)<br>
*>    A*x = lambda*B*x<br>
*> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the<br>
*> alternate form of the GNEP<br>
*>    mu*A*y = B*y.<br>
*> Real eigenvalues can be read directly from the generalized Schur<br>
*> form: <br>
*>   alpha = S(i,i), beta = P(i,i).<br>
*><br>
*> Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix<br>
*>      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),<br>
*>      pp. 241--256.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>          = 'E': Compute eigenvalues only;<br>
*>          = 'S': Compute eigenvalues and the Schur form. <br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPQ<br>
*> \verbatim<br>
*>          COMPQ is CHARACTER*1<br>
*>          = 'N': Left Schur vectors (Q) are not computed;<br>
*>          = 'I': Q is initialized to the unit matrix and the matrix Q<br>
*>                 of left Schur vectors of (H,T) is returned;<br>
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry and<br>
*>                 the product Q1*Q is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>          = 'N': Right Schur vectors (Z) are not computed;<br>
*>          = 'I': Z is initialized to the unit matrix and the matrix Z<br>
*>                 of right Schur vectors of (H,T) is returned;<br>
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry and<br>
*>                 the product Z1*Z is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrices H, T, Q, and Z.  N >= 0.<br>
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
*>          ILO and IHI mark the rows and columns of H which are in<br>
*>          Hessenberg form.  It is assumed that A is already upper<br>
*>          triangular in rows and columns 1:ILO-1 and IHI+1:N.<br>
*>          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is DOUBLE PRECISION array, dimension (LDH, N)<br>
*>          On entry, the N-by-N upper Hessenberg matrix H.<br>
*>          On exit, if JOB = 'S', H contains the upper quasi-triangular<br>
*>          matrix S from the generalized Schur factorization.<br>
*>          If JOB = 'E', the diagonal blocks of H match those of S, but<br>
*>          the rest of H is unspecified.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H.  LDH >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] T<br>
*> \verbatim<br>
*>          T is DOUBLE PRECISION array, dimension (LDT, N)<br>
*>          On entry, the N-by-N upper triangular matrix T.<br>
*>          On exit, if JOB = 'S', T contains the upper triangular<br>
*>          matrix P from the generalized Schur factorization;<br>
*>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S<br>
*>          are reduced to positive diagonal form, i.e., if H(j+1,j) is<br>
*>          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and<br>
*>          T(j+1,j+1) > 0.<br>
*>          If JOB = 'E', the diagonal blocks of T match those of P, but<br>
*>          the rest of T is unspecified.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDT<br>
*> \verbatim<br>
*>          LDT is INTEGER<br>
*>          The leading dimension of the array T.  LDT >= max( 1, N ).<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAR<br>
*> \verbatim<br>
*>          ALPHAR is DOUBLE PRECISION array, dimension (N)<br>
*>          The real parts of each scalar alpha defining an eigenvalue<br>
*>          of GNEP.<br>
*> \endverbatim<br>
*><br>
*> \param[out] ALPHAI<br>
*> \verbatim<br>
*>          ALPHAI is DOUBLE PRECISION array, dimension (N)<br>
*>          The imaginary parts of each scalar alpha defining an<br>
*>          eigenvalue of GNEP.<br>
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if<br>
*>          positive, then the j-th and (j+1)-st eigenvalues are a<br>
*>          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).<br>
*> \endverbatim<br>
*><br>
*> \param[out] BETA<br>
*> \verbatim<br>
*>          BETA is DOUBLE PRECISION array, dimension (N)<br>
*>          The scalars beta that define the eigenvalues of GNEP.<br>
*>          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and<br>
*>          beta = BETA(j) represent the j-th eigenvalue of the matrix<br>
*>          pair (A,B), in one of the forms lambda = alpha/beta or<br>
*>          mu = beta/alpha.  Since either lambda or mu may overflow,<br>
*>          they should not, in general, be computed.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Q<br>
*> \verbatim<br>
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)<br>
*>          On entry, if COMPQ = 'V', the orthogonal matrix Q1 used in<br>
*>          the reduction of (A,B) to generalized Hessenberg form.<br>
*>          On exit, if COMPQ = 'I', the orthogonal matrix of left Schur<br>
*>          vectors of (H,T), and if COMPQ = 'V', the orthogonal matrix<br>
*>          of left Schur vectors of (A,B).<br>
*>          Not referenced if COMPQ = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDQ<br>
*> \verbatim<br>
*>          LDQ is INTEGER<br>
*>          The leading dimension of the array Q.  LDQ >= 1.<br>
*>          If COMPQ='V' or 'I', then LDQ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)<br>
*>          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in<br>
*>          the reduction of (A,B) to generalized Hessenberg form.<br>
*>          On exit, if COMPZ = 'I', the orthogonal matrix of<br>
*>          right Schur vectors of (H,T), and if COMPZ = 'V', the<br>
*>          orthogonal matrix of right Schur vectors of (A,B).<br>
*>          Not referenced if COMPZ = 'N'.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>          The leading dimension of the array Z.  LDZ >= 1.<br>
*>          If COMPZ='V' or 'I', then LDZ >= N.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))<br>
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>          The dimension of the array WORK.  LWORK >= max(1,N).<br>
*><br>
*>          If LWORK = -1, then a workspace query is assumed; the routine<br>
*>          only calculates the optimal size of the WORK array, returns<br>
*>          this value as the first entry of the WORK array, and no error<br>
*>          message related to LWORK is issued by XERBLA.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0: successful exit<br>
*>          < 0: if INFO = -i, the i-th argument had an illegal value<br>
*>          = 1,...,N: the QZ iteration did not converge.  (H,T) is not<br>
*>                     in Schur form, but ALPHAR(i), ALPHAI(i), and<br>
*>                     BETA(i), i=INFO+1,...,N should be correct.<br>
*>          = N+1,...,2*N: the shift calculation failed.  (H,T) is not<br>
*>                     in Schur form, but ALPHAR(i), ALPHAI(i), and<br>
*>                     BETA(i), i=INFO-N+1,...,N should be correct.<br>
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
*> \ingroup doubleGEcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Iteration counters:<br>
*><br>
*>  JITER  -- counts iterations.<br>
*>  IITER  -- counts iterations run since ILAST was last<br>
*>            changed.  This is therefore reset only when a 1-by-1 or<br>
*>            2-by-2 block deflates off the bottom.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dhgeqz_(CHARACTER JOB,CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] T,INTEGER LDT,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);
/**
*> \brief \b DHSEIN<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DHSEIN + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhsein.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhsein.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhsein.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI,<br>
*                          VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL,<br>
*                          IFAILR, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       CHARACTER          EIGSRC, INITV, SIDE<br>
*       INTEGER            INFO, LDH, LDVL, LDVR, M, MM, N<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       LOGICAL            SELECT( * )<br>
*       INTEGER            IFAILL( * ), IFAILR( * )<br>
*       DOUBLE PRECISION   H( LDH, * ), VL( LDVL, * ), VR( LDVR, * ),<br>
*      $                   WI( * ), WORK( * ), WR( * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*> DHSEIN uses inverse iteration to find specified right and/or left<br>
*> eigenvectors of a real upper Hessenberg matrix H.<br>
*><br>
*> The right eigenvector x and the left eigenvector y of the matrix H<br>
*> corresponding to an eigenvalue w are defined by:<br>
*><br>
*>              H * x = w * x,     y**h * H = w * y**h<br>
*><br>
*> where y**h denotes the conjugate transpose of the vector y.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] SIDE<br>
*> \verbatim<br>
*>          SIDE is CHARACTER*1<br>
*>          = 'R': compute right eigenvectors only;<br>
*>          = 'L': compute left eigenvectors only;<br>
*>          = 'B': compute both right and left eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in] EIGSRC<br>
*> \verbatim<br>
*>          EIGSRC is CHARACTER*1<br>
*>          Specifies the source of eigenvalues supplied in (WR,WI):<br>
*>          = 'Q': the eigenvalues were found using DHSEQR; thus, if<br>
*>                 H has zero subdiagonal elements, and so is<br>
*>                 block-triangular, then the j-th eigenvalue can be<br>
*>                 assumed to be an eigenvalue of the block containing<br>
*>                 the j-th row/column.  This property allows DHSEIN to<br>
*>                 perform inverse iteration on just one diagonal block.<br>
*>          = 'N': no assumptions are made on the correspondence<br>
*>                 between eigenvalues and diagonal blocks.  In this<br>
*>                 case, DHSEIN must always perform inverse iteration<br>
*>                 using the whole matrix H.<br>
*> \endverbatim<br>
*><br>
*> \param[in] INITV<br>
*> \verbatim<br>
*>          INITV is CHARACTER*1<br>
*>          = 'N': no initial vectors are supplied;<br>
*>          = 'U': user-supplied initial vectors are stored in the arrays<br>
*>                 VL and/or VR.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] SELECT<br>
*> \verbatim<br>
*>          SELECT is LOGICAL array, dimension (N)<br>
*>          Specifies the eigenvectors to be computed. To select the<br>
*>          real eigenvector corresponding to a real eigenvalue WR(j),<br>
*>          SELECT(j) must be set to .TRUE.. To select the complex<br>
*>          eigenvector corresponding to a complex eigenvalue<br>
*>          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),<br>
*>          either SELECT(j) or SELECT(j+1) or both must be set to<br>
*>          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is<br>
*>          .FALSE..<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>          The order of the matrix H.  N >= 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in] H<br>
*> \verbatim<br>
*>          H is DOUBLE PRECISION array, dimension (LDH,N)<br>
*>          The upper Hessenberg matrix H.<br>
*>          If a NaN is detected in H, the routine will return with INFO=-6.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>          The leading dimension of the array H.  LDH >= max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] WR<br>
*> \verbatim<br>
*>          WR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[in] WI<br>
*> \verbatim<br>
*>          WI is DOUBLE PRECISION array, dimension (N)<br>
*><br>
*>          On entry, the real and imaginary parts of the eigenvalues of<br>
*>          H; a complex conjugate pair of eigenvalues must be stored in<br>
*>          consecutive elements of WR and WI.<br>
*>          On exit, WR may have been altered since close eigenvalues<br>
*>          are perturbed slightly in searching for independent<br>
*>          eigenvectors.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VL<br>
*> \verbatim<br>
*>          VL is DOUBLE PRECISION array, dimension (LDVL,MM)<br>
*>          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must<br>
*>          contain starting vectors for the inverse iteration for the<br>
*>          left eigenvectors; the starting vector for each eigenvector<br>
*>          must be in the same column(s) in which the eigenvector will<br>
*>          be stored.<br>
*>          On exit, if SIDE = 'L' or 'B', the left eigenvectors<br>
*>          specified by SELECT will be stored consecutively in the<br>
*>          columns of VL, in the same order as their eigenvalues. A<br>
*>          complex eigenvector corresponding to a complex eigenvalue is<br>
*>          stored in two consecutive columns, the first holding the real<br>
*>          part and the second the imaginary part.<br>
*>          If SIDE = 'R', VL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVL<br>
*> \verbatim<br>
*>          LDVL is INTEGER<br>
*>          The leading dimension of the array VL.<br>
*>          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] VR<br>
*> \verbatim<br>
*>          VR is DOUBLE PRECISION array, dimension (LDVR,MM)<br>
*>          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must<br>
*>          contain starting vectors for the inverse iteration for the<br>
*>          right eigenvectors; the starting vector for each eigenvector<br>
*>          must be in the same column(s) in which the eigenvector will<br>
*>          be stored.<br>
*>          On exit, if SIDE = 'R' or 'B', the right eigenvectors<br>
*>          specified by SELECT will be stored consecutively in the<br>
*>          columns of VR, in the same order as their eigenvalues. A<br>
*>          complex eigenvector corresponding to a complex eigenvalue is<br>
*>          stored in two consecutive columns, the first holding the real<br>
*>          part and the second the imaginary part.<br>
*>          If SIDE = 'L', VR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDVR<br>
*> \verbatim<br>
*>          LDVR is INTEGER<br>
*>          The leading dimension of the array VR.<br>
*>          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.<br>
*> \endverbatim<br>
*><br>
*> \param[in] MM<br>
*> \verbatim<br>
*>          MM is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR. MM >= M.<br>
*> \endverbatim<br>
*><br>
*> \param[out] M<br>
*> \verbatim<br>
*>          M is INTEGER<br>
*>          The number of columns in the arrays VL and/or VR required to<br>
*>          store the eigenvectors; each selected real eigenvector<br>
*>          occupies one column and each selected complex eigenvector<br>
*>          occupies two columns.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension ((N+2)*N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAILL<br>
*> \verbatim<br>
*>          IFAILL is INTEGER array, dimension (MM)<br>
*>          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left<br>
*>          eigenvector in the i-th column of VL (corresponding to the<br>
*>          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the<br>
*>          eigenvector converged satisfactorily. If the i-th and (i+1)th<br>
*>          columns of VL hold a complex eigenvector, then IFAILL(i) and<br>
*>          IFAILL(i+1) are set to the same value.<br>
*>          If SIDE = 'R', IFAILL is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] IFAILR<br>
*> \verbatim<br>
*>          IFAILR is INTEGER array, dimension (MM)<br>
*>          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right<br>
*>          eigenvector in the i-th column of VR (corresponding to the<br>
*>          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the<br>
*>          eigenvector converged satisfactorily. If the i-th and (i+1)th<br>
*>          columns of VR hold a complex eigenvector, then IFAILR(i) and<br>
*>          IFAILR(i+1) are set to the same value.<br>
*>          If SIDE = 'L', IFAILR is not referenced.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>          = 0:  successful exit<br>
*>          < 0:  if INFO = -i, the i-th argument had an illegal value<br>
*>          > 0:  if INFO = i, i is the number of eigenvectors which<br>
*>                failed to converge; see IFAILL and IFAILR for further<br>
*>                details.<br>
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
*> \date November 2013<br>
*<br>
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>  Each eigenvector is normalized so that the element of largest<br>
*>  magnitude has magnitude 1; here the magnitude of a complex number<br>
*>  (x,y) is taken to be |x|+|y|.<br>
*> \endverbatim<br>
*><br>
*  =====================================================================<br>
*/
	public void dhsein_(CHARACTER SIDE,CHARACTER EIGSRC,CHARACTER INITV,boolean[] SELECT,INTEGER N,double[] H,INTEGER LDH,double[] WR,double[] WI,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,int[] IFAILL,int[] IFAILR,INTEGER INFO);
/**
*> \brief \b DHSEQR<br>
*<br>
*  =========== DOCUMENTATION ===========<br>
*<br>
* Online html documentation available at <br>
*            http://www.netlib.org/lapack/explore-html/ <br>
*<br>
*> \htmlonly<br>
*> Download DHSEQR + dependencies <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhseqr.f"> <br>
*> [TGZ]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhseqr.f"> <br>
*> [ZIP]</a> <br>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhseqr.f"> <br>
*> [TXT]</a><br>
*> \endhtmlonly <br>
*<br>
*  Definition:<br>
*  ===========<br>
*<br>
*       SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,<br>
*                          LDZ, WORK, LWORK, INFO )<br>
* <br>
*       .. Scalar Arguments ..<br>
*       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N<br>
*       CHARACTER          COMPZ, JOB<br>
*       ..<br>
*       .. Array Arguments ..<br>
*       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),<br>
*      $                   Z( LDZ, * )<br>
*       ..<br>
*  <br>
*<br>
*> \par Purpose:<br>
*  =============<br>
*><br>
*> \verbatim<br>
*><br>
*>    DHSEQR computes the eigenvalues of a Hessenberg matrix H<br>
*>    and, optionally, the matrices T and Z from the Schur decomposition<br>
*>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the<br>
*>    Schur form), and Z is the orthogonal matrix of Schur vectors.<br>
*><br>
*>    Optionally Z may be postmultiplied into an input orthogonal<br>
*>    matrix Q so that this routine can give the Schur factorization<br>
*>    of a matrix A which has been reduced to the Hessenberg form H<br>
*>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.<br>
*> \endverbatim<br>
*<br>
*  Arguments:<br>
*  ==========<br>
*<br>
*> \param[in] JOB<br>
*> \verbatim<br>
*>          JOB is CHARACTER*1<br>
*>           = 'E':  compute eigenvalues only;<br>
*>           = 'S':  compute eigenvalues and the Schur form T.<br>
*> \endverbatim<br>
*><br>
*> \param[in] COMPZ<br>
*> \verbatim<br>
*>          COMPZ is CHARACTER*1<br>
*>           = 'N':  no Schur vectors are computed;<br>
*>           = 'I':  Z is initialized to the unit matrix and the matrix Z<br>
*>                   of Schur vectors of H is returned;<br>
*>           = 'V':  Z must contain an orthogonal matrix Q on entry, and<br>
*>                   the product Q*Z is returned.<br>
*> \endverbatim<br>
*><br>
*> \param[in] N<br>
*> \verbatim<br>
*>          N is INTEGER<br>
*>           The order of the matrix H.  N .GE. 0.<br>
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
*><br>
*>           It is assumed that H is already upper triangular in rows<br>
*>           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally<br>
*>           set by a previous call to DGEBAL, and then passed to ZGEHRD<br>
*>           when the matrix output by DGEBAL is reduced to Hessenberg<br>
*>           form. Otherwise ILO and IHI should be set to 1 and N<br>
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.<br>
*>           If N = 0, then ILO = 1 and IHI = 0.<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] H<br>
*> \verbatim<br>
*>          H is DOUBLE PRECISION array, dimension (LDH,N)<br>
*>           On entry, the upper Hessenberg matrix H.<br>
*>           On exit, if INFO = 0 and JOB = 'S', then H contains the<br>
*>           upper quasi-triangular matrix T from the Schur decomposition<br>
*>           (the Schur form); 2-by-2 diagonal blocks (corresponding to<br>
*>           complex conjugate pairs of eigenvalues) are returned in<br>
*>           standard form, with H(i,i) = H(i+1,i+1) and<br>
*>           H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the<br>
*>           contents of H are unspecified on exit.  (The output value of<br>
*>           H when INFO.GT.0 is given under the description of INFO<br>
*>           below.)<br>
*><br>
*>           Unlike earlier versions of DHSEQR, this subroutine may<br>
*>           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1<br>
*>           or j = IHI+1, IHI+2, ... N.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDH<br>
*> \verbatim<br>
*>          LDH is INTEGER<br>
*>           The leading dimension of the array H. LDH .GE. max(1,N).<br>
*> \endverbatim<br>
*><br>
*> \param[out] WR<br>
*> \verbatim<br>
*>          WR is DOUBLE PRECISION array, dimension (N)<br>
*> \endverbatim<br>
*><br>
*> \param[out] WI<br>
*> \verbatim<br>
*>          WI is DOUBLE PRECISION array, dimension (N)<br>
*><br>
*>           The real and imaginary parts, respectively, of the computed<br>
*>           eigenvalues. If two eigenvalues are computed as a complex<br>
*>           conjugate pair, they are stored in consecutive elements of<br>
*>           WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and<br>
*>           WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in<br>
*>           the same order as on the diagonal of the Schur form returned<br>
*>           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2<br>
*>           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and<br>
*>           WI(i+1) = -WI(i).<br>
*> \endverbatim<br>
*><br>
*> \param[in,out] Z<br>
*> \verbatim<br>
*>          Z is DOUBLE PRECISION array, dimension (LDZ,N)<br>
*>           If COMPZ = 'N', Z is not referenced.<br>
*>           If COMPZ = 'I', on entry Z need not be set and on exit,<br>
*>           if INFO = 0, Z contains the orthogonal matrix Z of the Schur<br>
*>           vectors of H.  If COMPZ = 'V', on entry Z must contain an<br>
*>           N-by-N matrix Q, which is assumed to be equal to the unit<br>
*>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,<br>
*>           if INFO = 0, Z contains Q*Z.<br>
*>           Normally Q is the orthogonal matrix generated by DORGHR<br>
*>           after the call to DGEHRD which formed the Hessenberg matrix<br>
*>           H. (The output value of Z when INFO.GT.0 is given under<br>
*>           the description of INFO below.)<br>
*> \endverbatim<br>
*><br>
*> \param[in] LDZ<br>
*> \verbatim<br>
*>          LDZ is INTEGER<br>
*>           The leading dimension of the array Z.  if COMPZ = 'I' or<br>
*>           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.<br>
*> \endverbatim<br>
*><br>
*> \param[out] WORK<br>
*> \verbatim<br>
*>          WORK is DOUBLE PRECISION array, dimension (LWORK)<br>
*>           On exit, if INFO = 0, WORK(1) returns an estimate of<br>
*>           the optimal value for LWORK.<br>
*> \endverbatim<br>
*><br>
*> \param[in] LWORK<br>
*> \verbatim<br>
*>          LWORK is INTEGER<br>
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)<br>
*>           is sufficient and delivers very good and sometimes<br>
*>           optimal performance.  However, LWORK as large as 11*N<br>
*>           may be required for optimal performance.  A workspace<br>
*>           query is recommended to determine the optimal workspace<br>
*>           size.<br>
*><br>
*>           If LWORK = -1, then DHSEQR does a workspace query.<br>
*>           In this case, DHSEQR checks the input parameters and<br>
*>           estimates the optimal workspace size for the given<br>
*>           values of N, ILO and IHI.  The estimate is returned<br>
*>           in WORK(1).  No error message related to LWORK is<br>
*>           issued by XERBLA.  Neither H nor Z are accessed.<br>
*> \endverbatim<br>
*><br>
*> \param[out] INFO<br>
*> \verbatim<br>
*>          INFO is INTEGER<br>
*>             =  0:  successful exit<br>
*>           .LT. 0:  if INFO = -i, the i-th argument had an illegal<br>
*>                    value<br>
*>           .GT. 0:  if INFO = i, DHSEQR failed to compute all of<br>
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR<br>
*>                and WI contain those eigenvalues which have been<br>
*>                successfully computed.  (Failures are rare.)<br>
*><br>
*>                If INFO .GT. 0 and JOB = 'E', then on exit, the<br>
*>                remaining unconverged eigenvalues are the eigen-<br>
*>                values of the upper Hessenberg matrix rows and<br>
*>                columns ILO through INFO of the final, output<br>
*>                value of H.<br>
*><br>
*>                If INFO .GT. 0 and JOB   = 'S', then on exit<br>
*><br>
*>           (*)  (initial value of H)*U  = U*(final value of H)<br>
*><br>
*>                where U is an orthogonal matrix.  The final<br>
*>                value of H is upper Hessenberg and quasi-triangular<br>
*>                in rows and columns INFO+1 through IHI.<br>
*><br>
*>                If INFO .GT. 0 and COMPZ = 'V', then on exit<br>
*><br>
*>                  (final value of Z)  =  (initial value of Z)*U<br>
*><br>
*>                where U is the orthogonal matrix in (*) (regard-<br>
*>                less of the value of JOB.)<br>
*><br>
*>                If INFO .GT. 0 and COMPZ = 'I', then on exit<br>
*>                      (final value of Z)  = U<br>
*>                where U is the orthogonal matrix in (*) (regard-<br>
*>                less of the value of JOB.)<br>
*><br>
*>                If INFO .GT. 0 and COMPZ = 'N', then Z is not<br>
*>                accessed.<br>
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
*> \ingroup doubleOTHERcomputational<br>
*<br>
*> \par Contributors:<br>
*  ==================<br>
*><br>
*>       Karen Braman and Ralph Byers, Department of Mathematics,<br>
*>       University of Kansas, USA<br>
*<br>
*> \par Further Details:<br>
*  =====================<br>
*><br>
*> \verbatim<br>
*><br>
*>             Default values supplied by<br>
*>             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).<br>
*>             It is suggested that these defaults be adjusted in order<br>
*>             to attain best performance in each particular<br>
*>             computational environment.<br>
*><br>
*>            ISPEC=12: The DLAHQR vs DLAQR0 crossover point.<br>
*>                      Default: 75. (Must be at least 11.)<br>
*><br>
*>            ISPEC=13: Recommended deflation window size.<br>
*>                      This depends on ILO, IHI and NS.  NS is the<br>
*>                      number of simultaneous shifts returned<br>
*>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)<br>
*>                      The default for (IHI-ILO+1).LE.500 is NS.<br>
*>                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.<br>
*><br>
*>            ISPEC=14: Nibble crossover point. (See IPARMQ for<br>
*>                      details.)  Default: 14% of deflation window<br>
*>                      size.<br>
*><br>
*>            ISPEC=15: Number of simultaneous shifts in a multishift<br>
*>                      QR iteration.<br>
*><br>
*>                      If IHI-ILO+1 is ...<br>
*><br>
*>                      greater than      ...but less    ... the<br>
*>                      or equal to ...      than        default is<br>
*><br>
*>                           1               30          NS =   2(+)<br>
*>                          30               60          NS =   4(+)<br>
*>                          60              150          NS =  10(+)<br>
*>                         150              590          NS =  **<br>
*>                         590             3000          NS =  64<br>
*>                        3000             6000          NS = 128<br>
*>                        6000             infinity      NS = 256<br>
*><br>
*>                  (+)  By default some or all matrices of this order<br>
*>                       are passed to the implicit double shift routine<br>
*>                       DLAHQR and this parameter is ignored.  See<br>
*>                       ISPEC=12 above and comments in IPARMQ for<br>
*>                       details.<br>
*><br>
*>                 (**)  The asterisks (**) indicate an ad-hoc<br>
*>                       function of N increasing from 10 to 64.<br>
*><br>
*>            ISPEC=16: Select structured matrix multiply.<br>
*>                      If the number of simultaneous shifts (specified<br>
*>                      by ISPEC=15) is less than 14, then the default<br>
*>                      for ISPEC=16 is 0.  Otherwise the default for<br>
*>                      ISPEC=16 is 2.<br>
*> \endverbatim<br>
*<br>
*> \par References:<br>
*  ================<br>
*><br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3<br>
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages<br>
*>       929--947, 2002.<br>
*> \n<br>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR<br>
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal<br>
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.<br>
*<br>
*  =====================================================================<br>
*/
	public void dhseqr_(CHARACTER JOB,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] WR,double[] WI,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);

}