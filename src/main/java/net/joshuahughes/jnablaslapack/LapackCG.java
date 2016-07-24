package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCG extends Library
{

	public static LapackCG instance = (LapackCG) Native.loadLibrary("liblapack",LapackCG.class);

	public void cgbbrd_(CHARACTER VECT,INTEGER M,INTEGER N,INTEGER NCC,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] D,float[] E,float[] Q,INTEGER LDQ,float[] PT,INTEGER LDPT,float[] C,INTEGER LDC,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgbcon_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgbequ_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void cgbequb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void cgbrfs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgbrfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgbsv_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void cgbsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgbsvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL CLA_GBRPVGRW);
	public void cgbtf2_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
	public void cgbtrf_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
	public void cgbtrs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void cgebak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,float[] SCALE,INTEGER M,float[] V,INTEGER LDV,INTEGER INFO);
	public void cgebal_(CHARACTER JOB,INTEGER N,float[] A,INTEGER LDA,INTEGER ILO,INTEGER IHI,float[] SCALE,INTEGER INFO);
	public void cgebd2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] WORK,INTEGER INFO);
	public void cgebrd_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgecon_(CHARACTER NORM,INTEGER N,float[] A,INTEGER LDA,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgeequ_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void cgeequb_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void cgees_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,INTEGER N,float[] A,INTEGER LDA,INTEGER SDIM,float[] W,float[] VS,INTEGER LDVS,float[] WORK,INTEGER LWORK,float[] RWORK,boolean BWORK,INTEGER INFO);
	public void cgeesx_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,INTEGER SDIM,float[] W,float[] VS,INTEGER LDVS,REAL RCONDE,REAL RCONDV,float[] WORK,INTEGER LWORK,float[] RWORK,boolean BWORK,INTEGER INFO);
	public void cgeev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,float[] A,INTEGER LDA,float[] W,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cgeevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,float[] W,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,float[] SCALE,REAL ABNRM,float[] RCONDE,float[] RCONDV,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cgehd2_(INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void cgehrd_(INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgejsv_(CHARACTER JOBA,CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBR,CHARACTER JOBT,CHARACTER JOBP,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] SVA,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] CWORK,INTEGER LWORK,float[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER INFO);
	public void cgelq2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void cgelqf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgels_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgelsd_(INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] S,REAL RCOND,INTEGER RANK,float[] WORK,INTEGER LWORK,float[] RWORK,int[] IWORK,INTEGER INFO);
	public void cgelss_(INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] S,REAL RCOND,INTEGER RANK,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cgelsy_(INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,int[] JPVT,REAL RCOND,INTEGER RANK,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cgemqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER NB,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
	public void cgeql2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void cgeqlf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgeqp3_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cgeqr2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void cgeqr2p_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void cgeqrf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgeqrfp_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgeqrt_(INTEGER M,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] T,INTEGER LDT,float[] WORK,INTEGER INFO);
	public void cgeqrt2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] T,INTEGER LDT,INTEGER INFO);
	public void cgeqrt3_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] T,INTEGER LDT,INTEGER INFO);
	public void cgerfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgerfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgerq2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void cgerqf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgesc2_(INTEGER N,float[] A,INTEGER LDA,float[] RHS,int[] IPIV,int[] JPIV,REAL SCALE);
	public void cgesdd_(CHARACTER JOBZ,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] S,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] WORK,INTEGER LWORK,float[] RWORK,int[] IWORK,INTEGER INFO);
	public void cgesv_(INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void cgesvd_(CHARACTER JOBU,CHARACTER JOBVT,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] S,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cgesvdx_(CHARACTER JOBU,CHARACTER JOBVT,CHARACTER RANGE,INTEGER M,INTEGER N,float[] A,INTEGER LDA,REAL VL,REAL VU,INTEGER IL,INTEGER IU,INTEGER NS,float[] S,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] WORK,INTEGER LWORK,float[] RWORK,int[] IWORK,INTEGER INFO);
	public void cgesvj_(CHARACTER JOBA,CHARACTER JOBU,CHARACTER JOBV,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] SVA,INTEGER MV,float[] V,INTEGER LDV,float[] CWORK,INTEGER LWORK,float[] RWORK,INTEGER LRWORK,INTEGER INFO);
	public void cgesvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgesvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL CLA_GERPVGRW);
	public void cgetc2_(INTEGER N,float[] A,INTEGER LDA,int[] IPIV,int[] JPIV,INTEGER INFO);
	public void cgetf2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void cgetrf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void cgetrf2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void cgetri_(INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgetrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void cggbak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,float[] LSCALE,float[] RSCALE,INTEGER M,float[] V,INTEGER LDV,INTEGER INFO);
	public void cggbal_(CHARACTER JOB,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER ILO,INTEGER IHI,float[] LSCALE,float[] RSCALE,float[] WORK,INTEGER INFO);
	public void cgges_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER SDIM,float[] ALPHA,float[] BETA,float[] VSL,INTEGER LDVSL,float[] VSR,INTEGER LDVSR,float[] WORK,INTEGER LWORK,float[] RWORK,boolean BWORK,INTEGER INFO);
	public void cgges3_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER SDIM,float[] ALPHA,float[] BETA,float[] VSL,INTEGER LDVSL,float[] VSR,INTEGER LDVSR,float[] WORK,INTEGER LWORK,float[] RWORK,boolean BWORK,INTEGER INFO);
	public void cggesx_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER SDIM,float[] ALPHA,float[] BETA,float[] VSL,INTEGER LDVSL,float[] VSR,INTEGER LDVSR,float[] RCONDE,float[] RCONDV,float[] WORK,INTEGER LWORK,float[] RWORK,int[] IWORK,INTEGER LIWORK,boolean BWORK,INTEGER INFO);
	public void cggev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHA,float[] BETA,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cggev3_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHA,float[] BETA,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void cggevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHA,float[] BETA,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,float[] LSCALE,float[] RSCALE,REAL ABNRM,REAL BBNRM,float[] RCONDE,float[] RCONDV,float[] WORK,INTEGER LWORK,float[] RWORK,int[] IWORK,boolean BWORK,INTEGER INFO);
	public void cggglm_(INTEGER N,INTEGER M,INTEGER P,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] D,float[] X,float[] Y,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgghd3_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgghrd_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,INTEGER INFO);
	public void cgglse_(INTEGER M,INTEGER N,INTEGER P,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,float[] D,float[] X,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cggqrf_(INTEGER N,INTEGER M,INTEGER P,float[] A,INTEGER LDA,float[] TAUA,float[] B,INTEGER LDB,float[] TAUB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cggrqf_(INTEGER M,INTEGER P,INTEGER N,float[] A,INTEGER LDA,float[] TAUA,float[] B,INTEGER LDB,float[] TAUB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cggsvd3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER N,INTEGER P,INTEGER K,INTEGER L,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHA,float[] BETA,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] Q,INTEGER LDQ,float[] WORK,INTEGER LWORK,float[] RWORK,int[] IWORK,INTEGER INFO);
	public void cggsvp3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL TOLA,REAL TOLB,INTEGER K,INTEGER L,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] Q,INTEGER LDQ,int[] IWORK,float[] RWORK,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgsvj0_(CHARACTER JOBV,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] SVA,INTEGER MV,float[] V,INTEGER LDV,REAL EPS,REAL SFMIN,REAL TOL,INTEGER NSWEEP,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgsvj1_(CHARACTER JOBV,INTEGER M,INTEGER N,INTEGER N1,float[] A,INTEGER LDA,float[] D,float[] SVA,INTEGER MV,float[] V,INTEGER LDV,REAL EPS,REAL SFMIN,REAL TOL,INTEGER NSWEEP,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void cgtcon_(CHARACTER NORM,INTEGER N,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,INTEGER INFO);
	public void cgtrfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DLF,float[] DF,float[] DUF,float[] DU2,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgtsv_(INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] B,INTEGER LDB,INTEGER INFO);
	public void cgtsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DLF,float[] DF,float[] DUF,float[] DU2,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cgttrf_(INTEGER N,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,INTEGER INFO);
	public void cgttrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void cgtts2_(INTEGER ITRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,float[] B,INTEGER LDB);

}