package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDL extends Library
{

	public static LapackDL instance = (LapackDL) Native.loadLibrary("liblapack",LapackDL.class);

	public void dlabad_(DOUBLE SMALL,DOUBLE LARGE);
	public void dlabrd_(INTEGER M,INTEGER N,INTEGER NB,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAUQ,double[] TAUP,double[] X,INTEGER LDX,double[] Y,INTEGER LDY);
	public void dlacn2_(INTEGER N,double[] V,double[] X,int[] ISGN,DOUBLE EST,INTEGER KASE,int[] ISAVE);
	public void dlacon_(INTEGER N,double[] V,double[] X,int[] ISGN,DOUBLE EST,INTEGER KASE);
	public void dlacpy_(CHARACTER UPLO,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB);
	public void dladiv_(DOUBLE A,DOUBLE B,DOUBLE C,DOUBLE D,DOUBLE P,DOUBLE Q);
	public void dladiv1_(DOUBLE A,DOUBLE B,DOUBLE C,DOUBLE D,DOUBLE P,DOUBLE Q);
	public double dladiv2_(DOUBLE A,DOUBLE B,DOUBLE C,DOUBLE D,DOUBLE R,DOUBLE T);
	public void dlae2_(DOUBLE A,DOUBLE B,DOUBLE C,DOUBLE RT1,DOUBLE RT2);
	public void dlaebz_(INTEGER IJOB,INTEGER NITMAX,INTEGER N,INTEGER MMAX,INTEGER MINP,INTEGER NBMIN,DOUBLE ABSTOL,DOUBLE RELTOL,DOUBLE PIVMIN,double[] D,double[] E,double[] E2,int[] NVAL,double[] AB,double[] C,INTEGER MOUT,int[] NAB,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlaed0_(INTEGER ICOMPQ,INTEGER QSIZ,INTEGER N,double[] D,double[] E,double[] Q,INTEGER LDQ,double[] QSTORE,INTEGER LDQS,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlaed1_(INTEGER N,double[] D,double[] Q,INTEGER LDQ,int[] INDXQ,DOUBLE RHO,INTEGER CUTPNT,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlaed2_(INTEGER K,INTEGER N,INTEGER N1,double[] D,double[] Q,INTEGER LDQ,int[] INDXQ,DOUBLE RHO,double[] Z,double[] DLAMDA,double[] W,double[] Q2,int[] INDX,int[] INDXC,int[] INDXP,int[] COLTYP,INTEGER INFO);
	public void dlaed3_(INTEGER K,INTEGER N,INTEGER N1,double[] D,double[] Q,INTEGER LDQ,DOUBLE RHO,double[] DLAMDA,double[] Q2,int[] INDX,int[] CTOT,double[] W,double[] S,INTEGER INFO);
	public void dlaed4_(INTEGER N,INTEGER I,double[] D,double[] Z,double[] DELTA,DOUBLE RHO,DOUBLE DLAM,INTEGER INFO);
	public void dlaed5_(INTEGER I,double[] D,double[] Z,double[] DELTA,DOUBLE RHO,DOUBLE DLAM);
	public void dlaed6_(INTEGER KNITER,LOGICAL ORGATI,DOUBLE RHO,double[] D,double[] Z,DOUBLE FINIT,DOUBLE TAU,INTEGER INFO);
	public void dlaed7_(INTEGER ICOMPQ,INTEGER N,INTEGER QSIZ,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,double[] D,double[] Q,INTEGER LDQ,int[] INDXQ,DOUBLE RHO,INTEGER CUTPNT,double[] QSTORE,int[] QPTR,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,double[] GIVNUM,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlaed8_(INTEGER ICOMPQ,INTEGER K,INTEGER N,INTEGER QSIZ,double[] D,double[] Q,INTEGER LDQ,int[] INDXQ,DOUBLE RHO,INTEGER CUTPNT,double[] Z,double[] DLAMDA,double[] Q2,INTEGER LDQ2,double[] W,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,double[] GIVNUM,int[] INDXP,int[] INDX,INTEGER INFO);
	public void dlaed9_(INTEGER K,INTEGER KSTART,INTEGER KSTOP,INTEGER N,double[] D,double[] Q,INTEGER LDQ,DOUBLE RHO,double[] DLAMDA,double[] W,double[] S,INTEGER LDS,INTEGER INFO);
	public void dlaeda_(INTEGER N,INTEGER TLVLS,INTEGER CURLVL,INTEGER CURPBM,int[] PRMPTR,int[] PERM,int[] GIVPTR,int[] GIVCOL,double[] GIVNUM,double[] Q,int[] QPTR,double[] Z,double[] ZTEMP,INTEGER INFO);
	public void dlaein_(LOGICAL RIGHTV,LOGICAL NOINIT,INTEGER N,double[] H,INTEGER LDH,DOUBLE WR,DOUBLE WI,double[] VR,double[] VI,double[] B,INTEGER LDB,double[] WORK,DOUBLE EPS3,DOUBLE SMLNUM,DOUBLE BIGNUM,INTEGER INFO);
	public void dlaev2_(DOUBLE A,DOUBLE B,DOUBLE C,DOUBLE RT1,DOUBLE RT2,DOUBLE CS1,DOUBLE SN1);
	public void dlaexc_(LOGICAL WANTQ,INTEGER N,double[] T,INTEGER LDT,double[] Q,INTEGER LDQ,INTEGER J1,INTEGER N1,INTEGER N2,double[] WORK,INTEGER INFO);
	public void dlag2_(double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE SAFMIN,DOUBLE SCALE1,DOUBLE SCALE2,DOUBLE WR1,DOUBLE WR2,DOUBLE WI);
	public void dlag2s_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,float[] SA,INTEGER LDSA,INTEGER INFO);
	public void dlags2_(LOGICAL UPPER,DOUBLE A1,DOUBLE A2,DOUBLE A3,DOUBLE B1,DOUBLE B2,DOUBLE B3,DOUBLE CSU,DOUBLE SNU,DOUBLE CSV,DOUBLE SNV,DOUBLE CSQ,DOUBLE SNQ);
	public void dlagtf_(INTEGER N,double[] A,DOUBLE LAMBDA,double[] B,double[] C,DOUBLE TOL,double[] D,int[] IN,INTEGER INFO);
	public void dlagtm_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,DOUBLE ALPHA,double[] DL,double[] D,double[] DU,double[] X,INTEGER LDX,DOUBLE BETA,double[] B,INTEGER LDB);
	public void dlagts_(INTEGER JOB,INTEGER N,double[] A,double[] B,double[] C,double[] D,int[] IN,double[] Y,DOUBLE TOL,INTEGER INFO);
	public void dlagv2_(double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,DOUBLE CSL,DOUBLE SNL,DOUBLE CSR,DOUBLE SNR);
	public void dlahqr_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] WR,double[] WI,INTEGER ILOZ,INTEGER IHIZ,double[] Z,INTEGER LDZ,INTEGER INFO);
	public void dlahr2_(INTEGER N,INTEGER K,INTEGER NB,double[] A,INTEGER LDA,double[] TAU,double[] T,INTEGER LDT,double[] Y,INTEGER LDY);
	public void dlaic1_(INTEGER JOB,INTEGER J,double[] X,DOUBLE SEST,double[] W,DOUBLE GAMMA,DOUBLE SESTPR,DOUBLE S,DOUBLE C);
	public boolean dlaisnan_(DOUBLE DIN1,DOUBLE DIN2);
	public void dlaln2_(LOGICAL LTRANS,INTEGER NA,INTEGER NW,DOUBLE SMIN,DOUBLE CA,double[] A,INTEGER LDA,DOUBLE D1,DOUBLE D2,double[] B,INTEGER LDB,DOUBLE WR,DOUBLE WI,double[] X,INTEGER LDX,DOUBLE SCALE,DOUBLE XNORM,INTEGER INFO);
	public void dlals0_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER NRHS,double[] B,INTEGER LDB,double[] BX,INTEGER LDBX,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,double[] GIVNUM,INTEGER LDGNUM,double[] POLES,double[] DIFL,double[] DIFR,double[] Z,INTEGER K,DOUBLE C,DOUBLE S,double[] WORK,INTEGER INFO);
	public void dlalsa_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,double[] B,INTEGER LDB,double[] BX,INTEGER LDBX,double[] U,INTEGER LDU,double[] VT,int[] K,double[] DIFL,double[] DIFR,double[] Z,double[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,double[] GIVNUM,double[] C,double[] S,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlalsd_(CHARACTER UPLO,INTEGER SMLSIZ,INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB,DOUBLE RCOND,INTEGER RANK,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlamrg_(INTEGER N1,INTEGER N2,double[] A,INTEGER DTRD1,INTEGER DTRD2,int[] INDEX);
	public int dlaneg_(INTEGER N,double[] D,double[] LLD,DOUBLE SIGMA,DOUBLE PIVMIN,INTEGER R);
	public double dlangb_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] WORK);
	public double dlange_(CHARACTER NORM,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] WORK);
	public double dlangt_(CHARACTER NORM,INTEGER N,double[] DL,double[] D,double[] DU);
	public double dlanhs_(CHARACTER NORM,INTEGER N,double[] A,INTEGER LDA,double[] WORK);
	public double dlansb_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,INTEGER K,double[] AB,INTEGER LDAB,double[] WORK);
	public double dlansf_(CHARACTER NORM,CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,double[] WORK);
	public double dlansp_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,double[] AP,double[] WORK);
	public double dlanst_(CHARACTER NORM,INTEGER N,double[] D,double[] E);
	public double dlansy_(CHARACTER NORM,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] WORK);
	public double dlantb_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER K,double[] AB,INTEGER LDAB,double[] WORK);
	public double dlantp_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] AP,double[] WORK);
	public double dlantr_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] WORK);
	public void dlanv2_(DOUBLE A,DOUBLE B,DOUBLE C,DOUBLE D,DOUBLE RT1R,DOUBLE RT1I,DOUBLE RT2R,DOUBLE RT2I,DOUBLE CS,DOUBLE SN);
	public void dlapll_(INTEGER N,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,DOUBLE SSMIN);
	public void dlapmr_(LOGICAL FORWRD,INTEGER M,INTEGER N,double[] X,INTEGER LDX,int[] K);
	public void dlapmt_(LOGICAL FORWRD,INTEGER M,INTEGER N,double[] X,INTEGER LDX,int[] K);
	public double dlapy2_(DOUBLE X,DOUBLE Y);
	public double dlapy3_(DOUBLE X,DOUBLE Y,DOUBLE Z);
	public void dlaqgb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,CHARACTER EQUED);
	public void dlaqge_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,CHARACTER EQUED);
	public void dlaqp2_(INTEGER M,INTEGER N,INTEGER OFFSET,double[] A,INTEGER LDA,int[] JPVT,double[] TAU,double[] VN1,double[] VN2,double[] WORK);
	public void dlaqps_(INTEGER M,INTEGER N,INTEGER OFFSET,INTEGER NB,INTEGER KB,double[] A,INTEGER LDA,int[] JPVT,double[] TAU,double[] VN1,double[] VN2,double[] AUXV,double[] F,INTEGER LDF);
	public void dlaqr0_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] WR,double[] WI,INTEGER ILOZ,INTEGER IHIZ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dlaqr1_(INTEGER N,double[] H,INTEGER LDH,DOUBLE SR1,DOUBLE SI1,DOUBLE SR2,DOUBLE SI2,double[] V);
	public void dlaqr2_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,double[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,double[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,double[] SR,double[] SI,double[] V,INTEGER LDV,INTEGER NH,double[] T,INTEGER LDT,INTEGER NV,double[] WV,INTEGER LDWV,double[] WORK,INTEGER LWORK);
	public void dlaqr3_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NW,double[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,double[] Z,INTEGER LDZ,INTEGER NS,INTEGER ND,double[] SR,double[] SI,double[] V,INTEGER LDV,INTEGER NH,double[] T,INTEGER LDT,INTEGER NV,double[] WV,INTEGER LDWV,double[] WORK,INTEGER LWORK);
	public void dlaqr4_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] WR,double[] WI,INTEGER ILOZ,INTEGER IHIZ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dlaqr5_(LOGICAL WANTT,LOGICAL WANTZ,INTEGER KACC22,INTEGER N,INTEGER KTOP,INTEGER KBOT,INTEGER NSHFTS,double[] SR,double[] SI,double[] H,INTEGER LDH,INTEGER ILOZ,INTEGER IHIZ,double[] Z,INTEGER LDZ,double[] V,INTEGER LDV,double[] U,INTEGER LDU,INTEGER NV,double[] WV,INTEGER LDWV,INTEGER NH,double[] WH,INTEGER LDWH);
	public void dlaqsb_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] S,DOUBLE SCOND,DOUBLE AMAX,CHARACTER EQUED);
	public void dlaqsp_(CHARACTER UPLO,INTEGER N,double[] AP,double[] S,DOUBLE SCOND,DOUBLE AMAX,CHARACTER EQUED);
	public void dlaqsy_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,CHARACTER EQUED);
	public void dlaqtr_(LOGICAL LTRAN,LOGICAL LREAL,INTEGER N,double[] T,INTEGER LDT,double[] B,DOUBLE W,DOUBLE SCALE,double[] X,double[] WORK,INTEGER INFO);
	public void dlar1v_(INTEGER N,INTEGER B1,INTEGER BN,DOUBLE LAMBDA,double[] D,double[] L,double[] LD,double[] LLD,DOUBLE PIVMIN,DOUBLE GAPTOL,double[] Z,LOGICAL WANTNC,INTEGER NEGCNT,DOUBLE ZTZ,DOUBLE MINGMA,INTEGER R,int[] ISUPPZ,DOUBLE NRMINV,DOUBLE RESID,DOUBLE RQCORR,double[] WORK);
	public void dlar2v_(INTEGER N,double[] X,double[] Y,double[] Z,INTEGER INCX,double[] C,double[] S,INTEGER INCC);
	public void dlarf_(CHARACTER SIDE,INTEGER M,INTEGER N,double[] V,INTEGER INCV,DOUBLE TAU,double[] C,INTEGER LDC,double[] WORK);
	public void dlarfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] C,INTEGER LDC,double[] WORK,INTEGER LDWORK);
	public void dlarfg_(INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,DOUBLE TAU);
	public void dlarfgp_(INTEGER N,DOUBLE ALPHA,double[] X,INTEGER INCX,DOUBLE TAU);
	public void dlarft_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,double[] V,INTEGER LDV,double[] TAU,double[] T,INTEGER LDT);
	public void dlarfx_(CHARACTER SIDE,INTEGER M,INTEGER N,double[] V,DOUBLE TAU,double[] C,INTEGER LDC,double[] WORK);
	public void dlargv_(INTEGER N,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] C,INTEGER INCC);
	public void dlarnv_(INTEGER IDIST,int[] ISEED,INTEGER N,double[] X);
	public void dlarra_(INTEGER N,double[] D,double[] E,double[] E2,DOUBLE SPLTOL,DOUBLE TNRM,INTEGER NSPLIT,int[] ISPLIT,INTEGER INFO);
	public void dlarrb_(INTEGER N,double[] D,double[] LLD,INTEGER IFIRST,INTEGER ILAST,DOUBLE RTOL1,DOUBLE RTOL2,INTEGER OFFSET,double[] W,double[] WGAP,double[] WERR,double[] WORK,int[] IWORK,DOUBLE PIVMIN,DOUBLE SPDIAM,INTEGER TWIST,INTEGER INFO);
	public void dlarrc_(CHARACTER JOBT,INTEGER N,DOUBLE VL,DOUBLE VU,double[] D,double[] E,DOUBLE PIVMIN,INTEGER EIGCNT,INTEGER LCNT,INTEGER RCNT,INTEGER INFO);
	public void dlarrd_(CHARACTER RANGE,CHARACTER ORDER,INTEGER N,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,double[] GERS,DOUBLE RELTOL,double[] D,double[] E,double[] E2,DOUBLE PIVMIN,INTEGER NSPLIT,int[] ISPLIT,INTEGER M,double[] W,double[] WERR,DOUBLE WL,DOUBLE WU,int[] IBLOCK,int[] INDEXW,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlarre_(CHARACTER RANGE,INTEGER N,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,double[] D,double[] E,double[] E2,DOUBLE RTOL1,DOUBLE RTOL2,DOUBLE SPLTOL,INTEGER NSPLIT,int[] ISPLIT,INTEGER M,double[] W,double[] WERR,double[] WGAP,int[] IBLOCK,int[] INDEXW,double[] GERS,DOUBLE PIVMIN,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlarrf_(INTEGER N,double[] D,double[] L,double[] LD,INTEGER CLSTRT,INTEGER CLEND,double[] W,double[] WGAP,double[] WERR,DOUBLE SPDIAM,DOUBLE CLGAPL,DOUBLE CLGAPR,DOUBLE PIVMIN,DOUBLE SIGMA,double[] DPLUS,double[] LPLUS,double[] WORK,INTEGER INFO);
	public void dlarrj_(INTEGER N,double[] D,double[] E2,INTEGER IFIRST,INTEGER ILAST,DOUBLE RTOL,INTEGER OFFSET,double[] W,double[] WERR,double[] WORK,int[] IWORK,DOUBLE PIVMIN,DOUBLE SPDIAM,INTEGER INFO);
	public void dlarrk_(INTEGER N,INTEGER IW,DOUBLE GL,DOUBLE GU,double[] D,double[] E2,DOUBLE PIVMIN,DOUBLE RELTOL,DOUBLE W,DOUBLE WERR,INTEGER INFO);
	public void dlarrr_(INTEGER N,double[] D,double[] E,INTEGER INFO);
	public void dlarrv_(INTEGER N,DOUBLE VL,DOUBLE VU,double[] D,double[] L,DOUBLE PIVMIN,int[] ISPLIT,INTEGER M,INTEGER DOL,INTEGER DOU,DOUBLE MINRGP,DOUBLE RTOL1,DOUBLE RTOL2,double[] W,double[] WERR,double[] WGAP,int[] IBLOCK,int[] INDEXW,double[] GERS,double[] Z,INTEGER LDZ,int[] ISUPPZ,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlarscl2_(INTEGER M,INTEGER N,double[] D,double[] X,INTEGER LDX);
	public void dlartg_(DOUBLE F,DOUBLE G,DOUBLE CS,DOUBLE SN,DOUBLE R);
	public void dlartgp_(DOUBLE F,DOUBLE G,DOUBLE CS,DOUBLE SN,DOUBLE R);
	public void dlartgs_(DOUBLE X,DOUBLE Y,DOUBLE SIGMA,DOUBLE CS,DOUBLE SN);
	public void dlartv_(INTEGER N,double[] X,INTEGER INCX,double[] Y,INTEGER INCY,double[] C,double[] S,INTEGER INCC);
	public void dlaruv_(int[] ISEED,INTEGER N,double[] X);
	public void dlarz_(CHARACTER SIDE,INTEGER M,INTEGER N,INTEGER L,double[] V,INTEGER INCV,DOUBLE TAU,double[] C,INTEGER LDC,double[] WORK);
	public void dlarzb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] C,INTEGER LDC,double[] WORK,INTEGER LDWORK);
	public void dlarzt_(CHARACTER DIRECT,CHARACTER STOREV,INTEGER N,INTEGER K,double[] V,INTEGER LDV,double[] TAU,double[] T,INTEGER LDT);
	public void dlas2_(DOUBLE F,DOUBLE G,DOUBLE H,DOUBLE SSMIN,DOUBLE SSMAX);
	public void dlascl_(CHARACTER TYPE,INTEGER KL,INTEGER KU,DOUBLE CFROM,DOUBLE CTO,INTEGER M,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dlascl2_(INTEGER M,INTEGER N,double[] D,double[] X,INTEGER LDX);
	public void dlasd0_(INTEGER N,INTEGER SQRE,double[] D,double[] E,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,INTEGER SMLSIZ,int[] IWORK,double[] WORK,INTEGER INFO);
	public void dlasd1_(INTEGER NL,INTEGER NR,INTEGER SQRE,double[] D,DOUBLE ALPHA,DOUBLE BETA,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,int[] IDXQ,int[] IWORK,double[] WORK,INTEGER INFO);
	public void dlasd2_(INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,double[] D,double[] Z,DOUBLE ALPHA,DOUBLE BETA,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] DSIGMA,double[] U2,INTEGER LDU2,double[] VT2,INTEGER LDVT2,int[] IDXP,int[] IDX,int[] IDXC,int[] IDXQ,int[] COLTYP,INTEGER INFO);
	public void dlasd3_(INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,double[] D,double[] Q,INTEGER LDQ,double[] DSIGMA,double[] U,INTEGER LDU,double[] U2,INTEGER LDU2,double[] VT,INTEGER LDVT,double[] VT2,INTEGER LDVT2,int[] IDXC,int[] CTOT,double[] Z,INTEGER INFO);
	public void dlasd4_(INTEGER N,INTEGER I,double[] D,double[] Z,double[] DELTA,DOUBLE RHO,DOUBLE SIGMA,double[] WORK,INTEGER INFO);
	public void dlasd5_(INTEGER I,double[] D,double[] Z,double[] DELTA,DOUBLE RHO,DOUBLE DSIGMA,double[] WORK);
	public void dlasd6_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,double[] D,double[] VF,double[] VL,DOUBLE ALPHA,DOUBLE BETA,int[] IDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,double[] GIVNUM,INTEGER LDGNUM,double[] POLES,double[] DIFL,double[] DIFR,double[] Z,INTEGER K,DOUBLE C,DOUBLE S,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlasd7_(INTEGER ICOMPQ,INTEGER NL,INTEGER NR,INTEGER SQRE,INTEGER K,double[] D,double[] Z,double[] ZW,double[] VF,double[] VFW,double[] VL,double[] VLW,DOUBLE ALPHA,DOUBLE BETA,double[] DSIGMA,int[] IDX,int[] IDXP,int[] IDXQ,int[] PERM,INTEGER GIVPTR,int[] GIVCOL,INTEGER LDGCOL,double[] GIVNUM,INTEGER LDGNUM,DOUBLE C,DOUBLE S,INTEGER INFO);
	public void dlasd8_(INTEGER ICOMPQ,INTEGER K,double[] D,double[] Z,double[] VF,double[] VL,double[] DIFL,double[] DIFR,INTEGER LDDIFR,double[] DSIGMA,double[] WORK,INTEGER INFO);
	public void dlasda_(INTEGER ICOMPQ,INTEGER SMLSIZ,INTEGER N,INTEGER SQRE,double[] D,double[] E,double[] U,INTEGER LDU,double[] VT,int[] K,double[] DIFL,double[] DIFR,double[] Z,double[] POLES,int[] GIVPTR,int[] GIVCOL,INTEGER LDGCOL,int[] PERM,double[] GIVNUM,double[] C,double[] S,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dlasdq_(CHARACTER UPLO,INTEGER SQRE,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,double[] D,double[] E,double[] VT,INTEGER LDVT,double[] U,INTEGER LDU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
	public void dlasdt_(INTEGER N,INTEGER LVL,INTEGER ND,int[] INODE,int[] NDIML,int[] NDIMR,INTEGER MSUB);
	public void dlaset_(CHARACTER UPLO,INTEGER M,INTEGER N,DOUBLE ALPHA,DOUBLE BETA,double[] A,INTEGER LDA);
	public void dlasq1_(INTEGER N,double[] D,double[] E,double[] WORK,INTEGER INFO);
	public void dlasq2_(INTEGER N,double[] Z,INTEGER INFO);
	public void dlasq3_(INTEGER I0,INTEGER N0,double[] Z,INTEGER PP,DOUBLE DMIN,DOUBLE SIGMA,DOUBLE DESIG,DOUBLE QMAX,INTEGER NFAIL,INTEGER ITER,INTEGER NDIV,LOGICAL IEEE,INTEGER TTYPE,DOUBLE DMIN1,DOUBLE DMIN2,DOUBLE DN,DOUBLE DN1,DOUBLE DN2,DOUBLE G,DOUBLE TAU);
	public void dlasq4_(INTEGER I0,INTEGER N0,double[] Z,INTEGER PP,INTEGER N0IN,DOUBLE DMIN,DOUBLE DMIN1,DOUBLE DMIN2,DOUBLE DN,DOUBLE DN1,DOUBLE DN2,DOUBLE TAU,INTEGER TTYPE,DOUBLE G);
	public void dlasq5_(INTEGER I0,INTEGER N0,double[] Z,INTEGER PP,DOUBLE TAU,DOUBLE SIGMA,DOUBLE DMIN,DOUBLE DMIN1,DOUBLE DMIN2,DOUBLE DN,DOUBLE DNM1,DOUBLE DNM2,LOGICAL IEEE,DOUBLE EPS);
	public void dlasq6_(INTEGER I0,INTEGER N0,double[] Z,INTEGER PP,DOUBLE DMIN,DOUBLE DMIN1,DOUBLE DMIN2,DOUBLE DN,DOUBLE DNM1,DOUBLE DNM2);
	public void dlasr_(CHARACTER SIDE,CHARACTER PIVOT,CHARACTER DIRECT,INTEGER M,INTEGER N,double[] C,double[] S,double[] A,INTEGER LDA);
	public void dlasrt_(CHARACTER ID,INTEGER N,double[] D,INTEGER INFO);
	public void dlassq_(INTEGER N,double[] X,INTEGER INCX,DOUBLE SCALE,DOUBLE SUMSQ);
	public void dlasv2_(DOUBLE F,DOUBLE G,DOUBLE H,DOUBLE SSMIN,DOUBLE SSMAX,DOUBLE SNR,DOUBLE CSR,DOUBLE SNL,DOUBLE CSL);
	public void dlaswp_(INTEGER N,double[] A,INTEGER LDA,INTEGER K1,INTEGER K2,int[] IPIV,INTEGER INCX);
	public void dlasy2_(LOGICAL LTRANL,LOGICAL LTRANR,INTEGER ISGN,INTEGER N1,INTEGER N2,double[] TL,INTEGER LDTL,double[] TR,INTEGER LDTR,double[] B,INTEGER LDB,DOUBLE SCALE,double[] X,INTEGER LDX,DOUBLE XNORM,INTEGER INFO);
	public void dlasyf_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,double[] A,INTEGER LDA,int[] IPIV,double[] W,INTEGER LDW,INTEGER INFO);
	public void dlasyf_rook_(CHARACTER UPLO,INTEGER N,INTEGER NB,INTEGER KB,double[] A,INTEGER LDA,int[] IPIV,double[] W,INTEGER LDW,INTEGER INFO);
	public void dlat2s_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,float[] SA,INTEGER LDSA,INTEGER INFO);
	public void dlatbs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] X,DOUBLE SCALE,double[] CNORM,INTEGER INFO);
	public void dlatdf_(INTEGER IJOB,INTEGER N,double[] Z,INTEGER LDZ,double[] RHS,DOUBLE RDSUM,DOUBLE RDSCAL,int[] IPIV,int[] JPIV);
	public void dlatps_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,double[] AP,double[] X,DOUBLE SCALE,double[] CNORM,INTEGER INFO);
	public void dlatrd_(CHARACTER UPLO,INTEGER N,INTEGER NB,double[] A,INTEGER LDA,double[] E,double[] TAU,double[] W,INTEGER LDW);
	public void dlatrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,CHARACTER NORMIN,INTEGER N,double[] A,INTEGER LDA,double[] X,DOUBLE SCALE,double[] CNORM,INTEGER INFO);
	public void dlatrz_(INTEGER M,INTEGER N,INTEGER L,double[] A,INTEGER LDA,double[] TAU,double[] WORK);
	public void dlauu2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dlauum_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dla_gbamv_(INTEGER TRANS,INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,DOUBLE ALPHA,double[] AB,INTEGER LDAB,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY,DOUBLE DLAMCH);
	public double dla_gbrcond_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,INTEGER CMODE,double[] C,INTEGER INFO,double[] WORK,int[] IWORK);
	public void dla_gbrfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,LOGICAL COLEQU,double[] C,double[] B,INTEGER LDB,double[] Y,INTEGER LDY,double[] BERR_OUT,INTEGER N_NORMS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,double[] RES,double[] AYB,double[] DY,double[] Y_TAIL,DOUBLE RCOND,INTEGER ITHRESH,DOUBLE RTHRESH,DOUBLE DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,DOUBLE DLAMCH,CHARACTER CHLA_TRANSTYPE);
	public double dla_gbrpvgrw_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NCOLS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB);
	public void dla_geamv_(INTEGER TRANS,INTEGER M,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY,DOUBLE DLAMCH);
	public double dla_gercond_(CHARACTER TRANS,INTEGER N,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,INTEGER CMODE,double[] C,INTEGER INFO,double[] WORK,int[] IWORK);
	public void dla_gerfsx_extended_(INTEGER PREC_TYPE,INTEGER TRANS_TYPE,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,double[] C,double[] B,INTEGER LDB,double[] Y,INTEGER LDY,double[] BERR_OUT,INTEGER N_NORMS,double[] ERRS_N,double[] ERRS_C,double[] RES,double[] AYB,double[] DY,double[] Y_TAIL,DOUBLE RCOND,INTEGER ITHRESH,DOUBLE RTHRESH,DOUBLE DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,DOUBLE DLAMCH,CHARACTER CHLA_TRANSTYPE);
	public double dla_gerpvgrw_(INTEGER N,INTEGER NCOLS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF);
	public void dla_lin_berr_(INTEGER N,INTEGER NZ,INTEGER NRHS,double[] RES,double[] AYB,double[] BERR);
	public double dla_porcond_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,INTEGER CMODE,double[] C,INTEGER INFO,double[] WORK,int[] IWORK);
	public void dla_porfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,LOGICAL COLEQU,double[] C,double[] B,INTEGER LDB,double[] Y,INTEGER LDY,double[] BERR_OUT,INTEGER N_NORMS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,double[] RES,double[] AYB,double[] DY,double[] Y_TAIL,DOUBLE RCOND,INTEGER ITHRESH,DOUBLE RTHRESH,DOUBLE DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
	public double dla_porpvgrw_(CHARACTER UPLO,INTEGER NCOLS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,double[] WORK);
	public void dla_syamv_(INTEGER UPLO,INTEGER N,DOUBLE ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,DOUBLE BETA,double[] Y,INTEGER INCY,DOUBLE DLAMCH);
	public double dla_syrcond_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,INTEGER CMODE,double[] C,INTEGER INFO,double[] WORK,int[] IWORK);
	public void dla_syrfsx_extended_(INTEGER PREC_TYPE,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,LOGICAL COLEQU,double[] C,double[] B,INTEGER LDB,double[] Y,INTEGER LDY,double[] BERR_OUT,INTEGER N_NORMS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,double[] RES,double[] AYB,double[] DY,double[] Y_TAIL,DOUBLE RCOND,INTEGER ITHRESH,DOUBLE RTHRESH,DOUBLE DZ_UB,LOGICAL IGNORE_CWISE,INTEGER INFO,INTEGER ILAUPLO);
	public double dla_syrpvgrw_(CHARACTER UPLO,INTEGER N,INTEGER INFO,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] WORK);
	public void dla_wwaddw_(INTEGER N,double[] X,double[] Y,double[] W);

}