package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDG extends Library
{

	public static LapackDG instance = (LapackDG) Native.loadLibrary("liblapack",LapackDG.class);

	public void dgbbrd_(CHARACTER VECT,INTEGER M,INTEGER N,INTEGER NCC,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] D,double[] E,double[] Q,INTEGER LDQ,double[] PT,INTEGER LDPT,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
	public void dgbcon_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgbequ_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
	public void dgbequb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
	public void dgbrfs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgbrfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgbsv_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void dgbsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgbsvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE DLA_GBRPVGRW);
	public void dgbtf2_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
	public void dgbtrf_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,double[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
	public void dgbtrs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,double[] AB,INTEGER LDAB,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void dgebak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,double[] SCALE,INTEGER M,double[] V,INTEGER LDV,INTEGER INFO);
	public void dgebal_(CHARACTER JOB,INTEGER N,double[] A,INTEGER LDA,INTEGER ILO,INTEGER IHI,double[] SCALE,INTEGER INFO);
	public void dgebd2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAUQ,double[] TAUP,double[] WORK,INTEGER INFO);
	public void dgebrd_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] E,double[] TAUQ,double[] TAUP,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgecon_(CHARACTER NORM,INTEGER N,double[] A,INTEGER LDA,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgeequ_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
	public void dgeequb_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] R,double[] C,DOUBLE ROWCND,DOUBLE COLCND,DOUBLE AMAX,INTEGER INFO);
	public void dgees_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,INTEGER N,double[] A,INTEGER LDA,INTEGER SDIM,double[] WR,double[] WI,double[] VS,INTEGER LDVS,double[] WORK,INTEGER LWORK,boolean BWORK,INTEGER INFO);
	public void dgeesx_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,INTEGER SDIM,double[] WR,double[] WI,double[] VS,INTEGER LDVS,DOUBLE RCONDE,DOUBLE RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,boolean BWORK,INTEGER INFO);
	public void dgeev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,double[] A,INTEGER LDA,double[] WR,double[] WI,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgeevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,double[] WR,double[] WI,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,double[] SCALE,DOUBLE ABNRM,double[] RCONDE,double[] RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void dgehd2_(INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
	public void dgehrd_(INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgejsv_(CHARACTER JOBA,CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBR,CHARACTER JOBT,CHARACTER JOBP,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] SVA,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void dgelq2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
	public void dgelqf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgels_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgelsd_(INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] S,DOUBLE RCOND,INTEGER RANK,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void dgelss_(INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] S,DOUBLE RCOND,INTEGER RANK,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgelsy_(INTEGER M,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,int[] JPVT,DOUBLE RCOND,INTEGER RANK,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgemqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER NB,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
	public void dgeql2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
	public void dgeqlf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgeqp3_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] JPVT,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgeqr2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
	public void dgeqr2p_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
	public void dgeqrf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgeqrfp_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgeqrt_(INTEGER M,INTEGER N,INTEGER NB,double[] A,INTEGER LDA,double[] T,INTEGER LDT,double[] WORK,INTEGER INFO);
	public void dgeqrt2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] T,INTEGER LDT,INTEGER INFO);
	public void dgeqrt3_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] T,INTEGER LDT,INTEGER INFO);
	public void dgerfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgerfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgerq2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER INFO);
	public void dgerqf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgesc2_(INTEGER N,double[] A,INTEGER LDA,double[] RHS,int[] IPIV,int[] JPIV,DOUBLE SCALE);
	public void dgesdd_(CHARACTER JOBZ,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] S,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void dgesv_(INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void dgesvd_(CHARACTER JOBU,CHARACTER JOBVT,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] S,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgesvdx_(CHARACTER JOBU,CHARACTER JOBVT,CHARACTER RANGE,INTEGER M,INTEGER N,double[] A,INTEGER LDA,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,INTEGER NS,double[] S,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void dgesvj_(CHARACTER JOBA,CHARACTER JOBU,CHARACTER JOBV,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] SVA,INTEGER MV,double[] V,INTEGER LDV,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgesvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgesvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] R,double[] C,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE DLA_GERPVGRW);
	public void dgetc2_(INTEGER N,double[] A,INTEGER LDA,int[] IPIV,int[] JPIV,INTEGER INFO);
	public void dgetf2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void dgetrf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void dgetrf2_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void dgetri_(INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgetrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void dggbak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,double[] LSCALE,double[] RSCALE,INTEGER M,double[] V,INTEGER LDV,INTEGER INFO);
	public void dggbal_(CHARACTER JOB,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER ILO,INTEGER IHI,double[] LSCALE,double[] RSCALE,double[] WORK,INTEGER INFO);
	public void dgges_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER SDIM,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VSL,INTEGER LDVSL,double[] VSR,INTEGER LDVSR,double[] WORK,INTEGER LWORK,boolean BWORK,INTEGER INFO);
	public void dgges3_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER SDIM,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VSL,INTEGER LDVSL,double[] VSR,INTEGER LDVSR,double[] WORK,INTEGER LWORK,boolean BWORK,INTEGER INFO);
	public void dggesx_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER SDIM,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VSL,INTEGER LDVSL,double[] VSR,INTEGER LDVSR,double[] RCONDE,double[] RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,boolean BWORK,INTEGER INFO);
	public void dggev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dggev3_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dggevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,double[] LSCALE,double[] RSCALE,DOUBLE ABNRM,DOUBLE BBNRM,double[] RCONDE,double[] RCONDV,double[] WORK,INTEGER LWORK,int[] IWORK,boolean BWORK,INTEGER INFO);
	public void dggglm_(INTEGER N,INTEGER M,INTEGER P,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] D,double[] X,double[] Y,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgghd3_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgghrd_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER INFO);
	public void dgglse_(INTEGER M,INTEGER N,INTEGER P,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,double[] D,double[] X,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dggqrf_(INTEGER N,INTEGER M,INTEGER P,double[] A,INTEGER LDA,double[] TAUA,double[] B,INTEGER LDB,double[] TAUB,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dggrqf_(INTEGER M,INTEGER P,INTEGER N,double[] A,INTEGER LDA,double[] TAUA,double[] B,INTEGER LDB,double[] TAUB,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dggsvd3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER N,INTEGER P,INTEGER K,INTEGER L,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHA,double[] BETA,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] Q,INTEGER LDQ,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void dggsvp3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE TOLA,DOUBLE TOLB,INTEGER K,INTEGER L,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] Q,INTEGER LDQ,int[] IWORK,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgsvj0_(CHARACTER JOBV,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] D,double[] SVA,INTEGER MV,double[] V,INTEGER LDV,DOUBLE EPS,DOUBLE SFMIN,DOUBLE TOL,INTEGER NSWEEP,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgsvj1_(CHARACTER JOBV,INTEGER M,INTEGER N,INTEGER N1,double[] A,INTEGER LDA,double[] D,double[] SVA,INTEGER MV,double[] V,INTEGER LDV,DOUBLE EPS,DOUBLE SFMIN,DOUBLE TOL,INTEGER NSWEEP,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dgtcon_(CHARACTER NORM,INTEGER N,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgtrfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DLF,double[] DF,double[] DUF,double[] DU2,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgtsv_(INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] B,INTEGER LDB,INTEGER INFO);
	public void dgtsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DLF,double[] DF,double[] DUF,double[] DU2,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dgttrf_(INTEGER N,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,INTEGER INFO);
	public void dgttrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void dgtts2_(INTEGER ITRANS,INTEGER N,INTEGER NRHS,double[] DL,double[] D,double[] DU,double[] DU2,int[] IPIV,double[] B,INTEGER LDB);

}