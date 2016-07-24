package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSL extends Library
{

	public static LapackSL instance = (LapackSL) Native.loadLibrary("liblapack",LapackSL.class);

	public void slabad_(REAL SMALL,REAL LARGE);
	public void slabrd_(INTEGER M,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] X,INTEGER LDX,float[] Y,INTEGER LDY);
	public void slacn2_(INTEGER N,float[] V,float[] X,int[] ISGN,REAL EST,INTEGER KASE,int[] ISAVE);
	public void slacon_(INTEGER N,float[] V,float[] X,int[] ISGN,REAL EST,INTEGER KASE);
	public void slacpy_(CHARACTER UPLO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB);
	public void sladiv_(REAL A,REAL B,REAL C,REAL D,REAL P,REAL Q);
	public void sladiv1_(REAL A,REAL B,REAL C,REAL D,REAL P,REAL Q);
	public float sladiv2_(REAL A,REAL B,REAL C,REAL D,REAL R,REAL T);
	public void slae2_(REAL A,REAL B,REAL C,REAL RT1,REAL RT2);
	public void slaebz_(INTEGER IJOB,INTEGER NITMAX,INTEGER N,INTEGER MMAX,INTEGER MINP,INTEGER NBMIN,REAL ABSTOL,REAL RELTOL,REAL PIVMIN,float[] D,float[] E,float[] E2,int[] NVAL,float[] AB,float[] C,INTEGER MOUT,int[] NAB,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slaed0_(INTEGER ICOMPQ,INTEGER QSIZ,INTEGER N,float[] D,float[] E,float[] Q,INTEGER LDQ,float[] QSTORE,INTEGER LDQS,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slaed1_(INTEGER N,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,INTEGER CUTPNT,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slaed2_(INTEGER K,INTEGER N,INTEGER N1,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,float[] Z,float[] DLAMDA,float[] W,float[] Q2,int[] INDX,int[] INDXC,int[] INDXP,int[] COLTYP,INTEGER INFO);
	public void slaed3_(INTEGER K,INTEGER N,INTEGER N1,float[] D,float[] Q,INTEGER LDQ,REAL RHO,float[] DLAMDA,float[] Q2,int[] INDX,int[] CTOT,float[] W,float[] S,INTEGER INFO);
	public void slaed4_(INTEGER N,INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL DLAM,INTEGER INFO);
	public void slaed5_(INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL DLAM);
	public void slaed6_(INTEGER KNITER,LOGICAL ORGATI,REAL RHO,float[] D,float[] Z,REAL FINIT,REAL TAU,INTEGER INFO);
	public void slaed7_(INTEGER ICOMPQ,INTEGER N,INTEGER QSIZ,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,INTEGER CUTPNT,float[] QSTORE,int[] QPTR,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,float[] GIVNUM,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slaed8_(INTEGER ICOMPQ,INTEGER K,INTEGER N,INTEGER QSIZ,float[] D,float[] Q,INTEGER LDQ,int[] INDXQ,REAL RHO,INTEGER CUTPNT,float[] Z,float[] DLAMDA,float[] Q2,INTEGER LDQ2,float[] W,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,float[] GIVNUM,int[] INDXP,int[] INDX,INTEGER INFO);
	public void slaed9_(INTEGER K,INTEGER KSTART,INTEGER KSTOP,INTEGER N,float[] D,float[] Q,INTEGER LDQ,REAL RHO,float[] DLAMDA,float[] W,float[] S,INTEGER LDS,INTEGER INFO);
	public void slaeda_(INTEGER N,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,float[] GIVNUM,float[] Q,int[] QPTR,float[] Z,float[] ZTEMP,INTEGER INFO);
	public void slaein_(LOGICAL RIGHTV,LOGICAL NOINIT,INTEGER N,float[] H,INTEGER LDH,REAL WR,REAL WI,float[] VR,float[] VI,float[] B,INTEGER LDB,float[] WORK,REAL EPS3,REAL SMLNUM,REAL BIGNUM,INTEGER INFO);
	public void slaev2_(REAL A,REAL B,REAL C,REAL RT1,REAL RT2,REAL CS1,REAL SN1);
	public void slaexc_(LOGICAL WANTQ,INTEGER N,float[] T,INTEGER LDT,float[] Q,INTEGER LDQ,INTEGER J1,INTEGER N1,INTEGER N2,float[] WORK,INTEGER INFO);
	public void slag2_(float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL SAFMIN,REAL SCALE1,REAL SCALE2,REAL WR1,REAL WR2,REAL WI);
	public void slag2d_(INTEGER M,INTEGER N,float[] SA,INTEGER LDSA,double[] A,INTEGER LDA,INTEGER INFO);
	public void slags2_(LOGICAL UPPER,REAL A1,REAL A2,REAL A3,REAL B1,REAL B2,REAL B3,REAL CSU,REAL SNU,REAL CSV,REAL SNV,REAL CSQ,REAL SNQ);
	public void slagtf_(INTEGER N,float[] A,REAL LAMBDA,float[] B,float[] C,REAL TOL,float[] D,int[] IN,INTEGER INFO);
	public void slagtm_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,REAL ALPHA,float[] DL,float[] D,float[] DU,float[] X,INTEGER LDX,REAL BETA,float[] B,INTEGER LDB);
	public void slagts_(INTEGER JOB,INTEGER N,float[] A,float[] B,float[] C,float[] D,int[] IN,float[] Y,REAL TOL,INTEGER INFO);
	public void slagv2_(float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHAR,float[] ALPHAI,float[] BETA,REAL CSL,REAL SNL,REAL CSR,REAL SNR);
	public void slahqr_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER INFO);
	public void slahr2_(INTEGER N,INTEGER K,INTEGER NB,float[] A,INTEGER LDA,float[] TAU,float[] T,INTEGER LDT,float[] Y,INTEGER LDY);
	public void slaic1_(INTEGER JOB,INTEGER J,float[] X,REAL SEST,float[] W,REAL GAMMA,REAL SESTPR,REAL S,REAL C);
	public boolean slaisnan_(REAL SIN1,REAL SIN2);
	public void slaln2_(LOGICAL LTRANS,INTEGER NA,INTEGER NW,REAL SMIN,REAL CA,float[] A,INTEGER LDA,REAL D1,REAL D2,float[] B,INTEGER LDB,REAL WR,REAL WI,float[] X,INTEGER LDX,REAL SCALE,REAL XNORM,INTEGER INFO);
	public void slals0_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER NRHS,float[] B,INTEGER LDB,float[] BX,INTEGER LDBX,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,float[] POLES,float[] DIFL,float[] DIFR,float[] Z,INTEGER K,REAL C,REAL S,float[] WORK,INTEGER INFO);
	public void slalsa_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,float[] B,INTEGER LDB,float[] BX,INTEGER LDBX,float[] U,INTEGER LDU,float[] VT,int[] K,float[] DIFL,float[] DIFR,float[] Z,float[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,float[] GIVNUM,float[] C,float[] S,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slalsd_(CHARACTER UPLO,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,REAL RCOND,INTEGER RANK,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slamrg_(INTEGER N1,INTEGER N2,float[] A,INTEGER STRD1,INTEGER STRD2,int[] INDEX);
	public int slaneg_(INTEGER N,float[] D,float[] LLD,REAL SIGMA,REAL PIVMIN,INTEGER R);
	public float slangb_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] WORK);
	public float slange_(CHARACTER NORM,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
	public float slangt_(CHARACTER NORM,INTEGER N,float[] DL,float[] D,float[] DU);
	public float slanhs_(CHARACTER NORM,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
	public float slansb_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
	public float slansf_(CHARACTER NORM,CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,float[] WORK);
	public float slansp_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] AP,float[] WORK);
	public float slanst_(CHARACTER NORM,INTEGER N,float[] D,float[] E);
	public float slansy_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
	public float slantb_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER K,float[] AB,INTEGER LDAB,float[] WORK);
	public float slantp_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,float[] AP,float[] WORK);
	public float slantr_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] WORK);
	public void slanv2_(REAL A,REAL B,REAL C,REAL D,REAL RT1R,REAL RT1I,REAL RT2R,REAL RT2I,REAL CS,REAL SN);
	public void slapll_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,REAL SSMIN);
	public void slapmr_(LOGICAL FORWRD,INTEGER M,INTEGER N,float[] X,INTEGER LDX,int[] K);
	public void slapmt_(LOGICAL FORWRD,INTEGER M,INTEGER N,float[] X,INTEGER LDX,int[] K);
	public float slapy2_(REAL X,REAL Y);
	public float slapy3_(REAL X,REAL Y,REAL Z);
	public void slaqgb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,CHARACTER EQUED);
	public void slaqge_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,CHARACTER EQUED);
	public void slaqp2_(INTEGER M,INTEGER N,INTEGER OFFSET,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] VN1,float[] VN2,float[] WORK);
	public void slaqps_(INTEGER M,INTEGER N,INTEGER OFFSET,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] VN1,float[] VN2,float[] AUXV,float[] F,INTEGER LDF);
	public void slaqr0_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void slaqr1_(INTEGER N,float[] H,INTEGER LDH,REAL SR1,REAL SI1,REAL SR2,REAL SI2,float[] V);
	public void slaqr2_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,float[] SR,float[] SI,float[] V,INTEGER LDV,INTEGER NH,float[] T,INTEGER LDT,INTEGER NV,float[] WV,INTEGER LDWV,float[] WORK,INTEGER LWORK);
	public void slaqr3_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,float[] SR,float[] SI,float[] V,INTEGER LDV,INTEGER NH,float[] T,INTEGER LDT,INTEGER NV,float[] WV,INTEGER LDWV,float[] WORK,INTEGER LWORK);
	public void slaqr4_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void slaqr5_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER KACC22,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NSHFTS,float[] SR,float[] SI,float[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,float[] Z,INTEGER LDZ,float[] V,INTEGER LDV,float[] U,INTEGER LDU,INTEGER NV,float[] WV,INTEGER LDWV,INTEGER NH,float[] WH,INTEGER LDWH);
	public void slaqsb_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
	public void slaqsp_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
	public void slaqsy_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,CHARACTER EQUED);
	public void slaqtr_(LOGICAL LTRAN,LOGICAL LREAL,INTEGER N,float[] T,INTEGER LDT,float[] B,REAL W,REAL SCALE,float[] X,float[] WORK,INTEGER INFO);
	public void slar1v_(INTEGER N,INTEGER B1,INTEGER BN,REAL LAMBDA,float[] D,float[] L,float[] LD,float[] LLD,REAL PIVMIN,REAL GAPTOL,float[] Z,LOGICAL WANTNC,INTEGER NEGCNT,REAL ZTZ,REAL MINGMA,INTEGER R,int[] ISUPPZ,REAL NRMINV,REAL RESID,REAL RQCORR,float[] WORK);
	public void slar2v_(INTEGER N,float[] X,float[] Y,float[] Z,INTEGER INCX,float[] C,float[] S,INTEGER INCC);
	public void slarf_(CHARACTER SIDE,INTEGER M,INTEGER N,float[] V,INTEGER INCV,REAL TAU,float[] C,INTEGER LDC,float[] WORK);
	public void slarfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER LDWORK);
	public void slarfg_(INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,REAL TAU);
	public void slarfgp_(INTEGER N,REAL ALPHA,float[] X,INTEGER INCX,REAL TAU);
	public void slarft_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] TAU,float[] T,INTEGER LDT);
	public void slarfx_(CHARACTER SIDE,INTEGER M,INTEGER N,float[] V,REAL TAU,float[] C,INTEGER LDC,float[] WORK);
	public void slargv_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] C,INTEGER INCC);
	public void slarnv_(INTEGER IDIST,int[] ISEED,INTEGER N,float[] X);
	public void slarra_(INTEGER N,float[] D,float[] E,float[] E2,REAL SPLTOL,REAL TNRM,INTEGER NSPLIT,int[] ISPLIT,INTEGER INFO);
	public void slarrb_(INTEGER N,float[] D,float[] LLD,INTEGER IFIRST,INTEGER ILAST,REAL RTOL1,REAL RTOL2,INTEGER OFFSET,float[] W,float[] WGAP,float[] WERR,float[] WORK,int[] IWORK,REAL PIVMIN,REAL SPDIAM,INTEGER TWIST,INTEGER INFO);
	public void slarrc_(CHARACTER JOBT,INTEGER N,REAL VL,REAL VU,float[] D,float[] E,REAL PIVMIN,INTEGER EIGCNT,INTEGER LCNT,INTEGER RCNT,INTEGER INFO);
	public void slarrd_(CHARACTER RANGE,CHARACTER ORDER,INTEGER N,REAL VL,REAL VU,INTEGER IL,INTEGER IU,float[] GERS,REAL RELTOL,float[] D,float[] E,float[] E2,REAL PIVMIN,INTEGER NSPLIT,int[] ISPLIT,INTEGER M,float[] W,float[] WERR,REAL WL,REAL WU,int[] IBLOCK,int[] INDEXW,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slarre_(CHARACTER RANGE,INTEGER N,REAL VL,REAL VU,INTEGER IL,INTEGER IU,float[] D,float[] E,float[] E2,REAL RTOL1,REAL RTOL2,REAL SPLTOL,INTEGER NSPLIT,int[] ISPLIT,INTEGER M,float[] W,float[] WERR,float[] WGAP,int[] IBLOCK,int[] INDEXW,float[] GERS,REAL PIVMIN,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slarrf_(INTEGER N,float[] D,float[] L,float[] LD,INTEGER CLSTRT,INTEGER CLEND,float[] W,float[] WGAP,float[] WERR,REAL SPDIAM,REAL CLGAPL,REAL CLGAPR,REAL PIVMIN,REAL SIGMA,float[] DPLUS,float[] LPLUS,float[] WORK,INTEGER INFO);
	public void slarrj_(INTEGER N,float[] D,float[] E2,INTEGER IFIRST,INTEGER ILAST,REAL RTOL,INTEGER OFFSET,float[] W,float[] WERR,float[] WORK,int[] IWORK,REAL PIVMIN,REAL SPDIAM,INTEGER INFO);
	public void slarrk_(INTEGER N,INTEGER IW,REAL GL,REAL GU,float[] D,float[] E2,REAL PIVMIN,REAL RELTOL,REAL W,REAL WERR,INTEGER INFO);
	public void slarrr_(INTEGER N,float[] D,float[] E,INTEGER INFO);
	public void slarrv_(INTEGER N,REAL VL,REAL VU,float[] D,float[] L,REAL PIVMIN,int[] ISPLIT,INTEGER M,INTEGER DOL,INTEGER DOU,REAL MINRGP,REAL RTOL1,REAL RTOL2,float[] W,float[] WERR,float[] WGAP,int[] IBLOCK,int[] INDEXW,float[] GERS,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slarscl2_(INTEGER M,INTEGER N,float[] D,float[] X,INTEGER LDX);
	public void slartg_(REAL F,REAL G,REAL CS,REAL SN,REAL R);
	public void slartgp_(REAL F,REAL G,REAL CS,REAL SN,REAL R);
	public void slartgs_(REAL X,REAL Y,REAL SIGMA,REAL CS,REAL SN);
	public void slartv_(INTEGER N,float[] X,INTEGER INCX,float[] Y,INTEGER INCY,float[] C,float[] S,INTEGER INCC);
	public void slaruv_(int[] ISEED,INTEGER N,float[] X);
	public void slarz_(CHARACTER SIDE,INTEGER M,INTEGER N,INTEGER L,float[] V,INTEGER INCV,REAL TAU,float[] C,INTEGER LDC,float[] WORK);
	public void slarzb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER LDWORK);
	public void slarzt_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,float[] V,INTEGER LDV,float[] TAU,float[] T,INTEGER LDT);
	public void slas2_(REAL F,REAL G,REAL H,REAL SSMIN,REAL SSMAX);
	public void slascl_(CHARACTER TYPE,INTEGER KL,INTEGER KU,REAL CFROM,REAL CTO,INTEGER M,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void slascl2_(INTEGER M,INTEGER N,float[] D,float[] X,INTEGER LDX);
	public void slasd0_(INTEGER N,INTEGER SQRE,float[] D,float[] E,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,INTEGER SMLSIZ,int[] IWORK,float[] WORK,INTEGER INFO);
	public void slasd1_(INTEGER NL,INTEGER NR,INTEGER SQRE,float[] D,REAL ALPHA,REAL BETA,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,int[] IDXQ,int[] IWORK,float[] WORK,INTEGER INFO);
	public void slasd2_(INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,float[] D,float[] Z,REAL ALPHA,REAL BETA,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] DSIGMA,float[] U2,INTEGER LDU2,float[] VT2,INTEGER LDVT2,int[] IDXP,int[] IDX,int[] IDXC,int[] IDXQ,int[] COLTYP,INTEGER INFO);
	public void slasd3_(INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,float[] D,float[] Q,INTEGER LDQ,float[] DSIGMA,float[] U,INTEGER LDU,float[] U2,INTEGER LDU2,float[] VT,INTEGER LDVT,float[] VT2,INTEGER LDVT2,int[] IDXC,int[] CTOT,float[] Z,INTEGER INFO);
	public void slasd4_(INTEGER N,INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL SIGMA,float[] WORK,INTEGER INFO);
	public void slasd5_(INTEGER I,float[] D,float[] Z,float[] DELTA,REAL RHO,REAL DSIGMA,float[] WORK);
	public void slasd6_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,float[] D,float[] VF,float[] VL,REAL ALPHA,REAL BETA,int[] IDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,float[] POLES,float[] DIFL,float[] DIFR,float[] Z,INTEGER K,REAL C,REAL S,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slasd7_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,float[] D,float[] Z,float[] ZW,float[] VF,float[] VFW,float[] VL,float[] VLW,REAL ALPHA,REAL BETA,float[] DSIGMA,int[] IDX,int[] IDXP,int[] IDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,float[] GIVNUM,INTEGER LDGNUM,REAL C,REAL S,INTEGER INFO);
	public void slasd8_(INTEGER ICOMPQ,INTEGER K,float[] D,float[] Z,float[] VF,float[] VL,float[] DIFL,float[] DIFR,INTEGER LDDIFR,float[] DSIGMA,float[] WORK,INTEGER INFO);
	public void slasda_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER SQRE,float[] D,float[] E,float[] U,INTEGER LDU,float[] VT,int[] K,float[] DIFL,float[] DIFR,float[] Z,float[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,float[] GIVNUM,float[] C,float[] S,float[] WORK,int[] IWORK,INTEGER INFO);
	public void slasdq_(CHARACTER UPLO,INTEGER SQRE,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,float[] D,float[] E,float[] VT,INTEGER LDVT,float[] U,INTEGER LDU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
	public void slasdt_(INTEGER N,INTEGER LVL,INTEGER ND,int[] INODE,int[] NDIML,int[] NDIMR,INTEGER MSUB);
	public void slaset_(CHARACTER UPLO,INTEGER M,INTEGER N,REAL ALPHA,REAL BETA,float[] A,INTEGER LDA);
	public void slasq1_(INTEGER N,float[] D,float[] E,float[] WORK,INTEGER INFO);
	public void slasq2_(INTEGER N,float[] Z,INTEGER INFO);
	public void slasq3_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,REAL DMIN,REAL SIGMA,REAL DESIG,REAL QMAX,INTEGER NFAIL,INTEGER ITER,INTEGER NDIV,LOGICAL IEEE,INTEGER TTYPE,REAL DMIN1,REAL DMIN2,REAL DN,REAL DN1,REAL DN2,REAL G,REAL TAU);
	public void slasq4_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,INTEGER N0IN,REAL DMIN,REAL DMIN1,REAL DMIN2,REAL DN,REAL DN1,REAL DN2,REAL TAU,INTEGER TTYPE,REAL G);
	public void slasq5_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,REAL TAU,REAL SIGMA,REAL DMIN,REAL DMIN1,REAL DMIN2,REAL DN,REAL DNM1,REAL DNM2,LOGICAL IEEE,REAL EPS);
	public void slasq6_(INTEGER I0,INTEGER N0,float[] Z,INTEGER PP,REAL DMIN,REAL DMIN1,REAL DMIN2,REAL DN,REAL DNM1,REAL DNM2);
	public void slasr_(CHARACTER SIDE,CHARACTER PIVOT,CHARACTER DIRECT,INTEGER M,INTEGER N,float[] C,float[] S,float[] A,INTEGER LDA);
	public void slasrt_(CHARACTER ID,INTEGER N,float[] D,INTEGER INFO);
	public void slassq_(INTEGER N,float[] X,INTEGER INCX,REAL SCALE,REAL SUMSQ);
	public void slasv2_(REAL F,REAL G,REAL H,REAL SSMIN,REAL SSMAX,REAL SNR,REAL CSR,REAL SNL,REAL CSL);
	public void slaswp_(INTEGER N,float[] A,INTEGER LDA,INTEGER K1,INTEGER K2,int[] IPIV,INTEGER INCX);
	public void slasy2_(LOGICAL LTRANL,LOGICAL LTRANR,INTEGER ISGN,INTEGER N1,INTEGER N2,float[] TL,INTEGER LDTL,float[] TR,INTEGER LDTR,float[] B,INTEGER LDB,REAL SCALE,float[] X,INTEGER LDX,REAL XNORM,INTEGER INFO);
	public void slasyf_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
	public void slasyf_rook_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,float[] A,INTEGER LDA,int[] IPIV,float[] W,INTEGER LDW,INTEGER INFO);
	public void slatbs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
	public void slatdf_(INTEGER IJOB,INTEGER N,float[] Z,INTEGER LDZ,float[] RHS,REAL RDSUM,REAL RDSCAL,int[] IPIV,int[] JPIV);
	public void slatps_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,float[] AP,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
	public void slatrd_(CHARACTER UPLO,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] E,float[] TAU,float[] W,INTEGER LDW);
	public void slatrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,float[] A,INTEGER LDA,float[] X,REAL SCALE,float[] CNORM,INTEGER INFO);
	public void slatrz_(INTEGER M,INTEGER N,INTEGER L,float[] A,INTEGER LDA,float[] TAU,float[] WORK);
	public void slauu2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void slauum_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void sla_gbamv_(INTEGER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,REAL ALPHA,float[] AB,INTEGER LDAB,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
	public float sla_gbrcond_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
	public void sla_gbrfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,REAL SLAMCH,CHARACTER CHLA_TRANSTYPE);
	public float sla_gbrpvgrw_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NCOLS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB);
	public void sla_geamv_(INTEGER TRANS,INTEGER M,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
	public float sla_gercond_(CHARACTER TRANS,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
	public void sla_gerfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERRS_N,float[] ERRS_C,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,REAL SLAMCH,CHARACTER CHLA_TRANSTYPE);
	public float sla_gerpvgrw_(INTEGER N,INTEGER NCOLS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF);
	public void sla_lin_berr_(INTEGER N,INTEGER NZ,INTEGER NRHS,float[] RES,float[] AYB,float[] BERR);
	public float sla_porcond_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
	public void sla_porfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
	public float sla_porpvgrw_(CHARACTER UPLO,INTEGER NCOLS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] WORK);
	public void sla_syamv_(INTEGER UPLO,INTEGER N,REAL ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,REAL BETA,float[] Y,INTEGER INCY,REAL SLAMCH);
	public float sla_syrcond_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,INTEGER CMODE,float[] C,INTEGER INFO,float[] WORK,int[] IWORK);
	public void sla_syrfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,float[] C,float[] B,INTEGER LDB,float[] Y,INTEGER LDY,float[] BERR_OUT,INTEGER N_NORMS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,float[] RES,float[] AYB,float[] DY,float[] Y_TAIL,REAL RCOND,INTEGER ITHRESH,REAL RTHRESH,REAL DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
	public float sla_syrpvgrw_(CHARACTER UPLO,INTEGER N,INTEGER INFO,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] WORK);
	public void sla_wwaddw_(INTEGER N,float[] X,float[] Y,float[] W);

}