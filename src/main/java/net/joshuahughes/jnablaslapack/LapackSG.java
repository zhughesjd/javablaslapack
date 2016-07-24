package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSG extends Library
{

	public static LapackSG instance = (LapackSG) Native.loadLibrary("liblapack",LapackSG.class);

	public void sgbbrd_(CHARACTER VECT,INTEGER M,INTEGER N,INTEGER NCC,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] D,float[] E,float[] Q,INTEGER LDQ,float[] PT,INTEGER LDPT,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
	public void sgbcon_(CHARACTER NORM,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgbequ_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void sgbequb_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void sgbrfs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgbrfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgbsv_(INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void sgbsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgbsvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL SLA_GBRPVGRW);
	public void sgbtf2_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
	public void sgbtrf_(INTEGER M,INTEGER N,INTEGER KL,INTEGER KU,float[] AB,INTEGER LDAB,int[] IPIV,INTEGER INFO);
	public void sgbtrs_(CHARACTER TRANS,INTEGER N,INTEGER KL,INTEGER KU,INTEGER NRHS,float[] AB,INTEGER LDAB,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void sgebak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,float[] SCALE,INTEGER M,float[] V,INTEGER LDV,INTEGER INFO);
	public void sgebal_(CHARACTER JOB,INTEGER N,float[] A,INTEGER LDA,INTEGER ILO,INTEGER IHI,float[] SCALE,INTEGER INFO);
	public void sgebd2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] WORK,INTEGER INFO);
	public void sgebrd_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] E,float[] TAUQ,float[] TAUP,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgecon_(CHARACTER NORM,INTEGER N,float[] A,INTEGER LDA,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgeequ_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void sgeequb_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] R,float[] C,REAL ROWCND,REAL COLCND,REAL AMAX,INTEGER INFO);
	public void sgees_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,INTEGER N,float[] A,INTEGER LDA,INTEGER SDIM,float[] WR,float[] WI,float[] VS,INTEGER LDVS,float[] WORK,INTEGER LWORK,boolean BWORK,INTEGER INFO);
	public void sgeesx_(CHARACTER JOBVS,CHARACTER SORT,LOGICAL SELECT,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,INTEGER SDIM,float[] WR,float[] WI,float[] VS,INTEGER LDVS,REAL RCONDE,REAL RCONDV,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,boolean BWORK,INTEGER INFO);
	public void sgeev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,float[] A,INTEGER LDA,float[] WR,float[] WI,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgeevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,float[] WR,float[] WI,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,float[] SCALE,REAL ABNRM,float[] RCONDE,float[] RCONDV,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void sgehd2_(INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void sgehrd_(INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgejsv_(CHARACTER JOBA,CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBR,CHARACTER JOBT,CHARACTER JOBP,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] SVA,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void sgelq2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void sgelqf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgels_(CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgelsd_(INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] S,REAL RCOND,INTEGER RANK,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void sgelss_(INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] S,REAL RCOND,INTEGER RANK,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgelsy_(INTEGER M,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,int[] JPVT,REAL RCOND,INTEGER RANK,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgemqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER NB,float[] V,INTEGER LDV,float[] T,INTEGER LDT,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
	public void sgeql2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void sgeqlf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgeqp3_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] JPVT,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgeqr2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void sgeqr2p_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void sgeqrf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgeqrfp_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgeqrt_(INTEGER M,INTEGER N,INTEGER NB,float[] A,INTEGER LDA,float[] T,INTEGER LDT,float[] WORK,INTEGER INFO);
	public void sgeqrt2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] T,INTEGER LDT,INTEGER INFO);
	public void sgeqrt3_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] T,INTEGER LDT,INTEGER INFO);
	public void sgerfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgerfsx_(CHARACTER TRANS,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgerq2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER INFO);
	public void sgerqf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgesc2_(INTEGER N,float[] A,INTEGER LDA,float[] RHS,int[] IPIV,int[] JPIV,REAL SCALE);
	public void sgesdd_(CHARACTER JOBZ,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] S,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void sgesv_(INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void sgesvd_(CHARACTER JOBU,CHARACTER JOBVT,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] S,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgesvdx_(CHARACTER JOBU,CHARACTER JOBVT,CHARACTER RANGE,INTEGER M,INTEGER N,float[] A,INTEGER LDA,REAL VL,REAL VU,INTEGER IL,INTEGER IU,INTEGER NS,float[] S,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void sgesvj_(CHARACTER JOBA,CHARACTER JOBU,CHARACTER JOBV,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] SVA,INTEGER MV,float[] V,INTEGER LDV,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgesvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgesvxx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,float[] R,float[] C,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL SLA_GERPVGRW);
	public void sgetc2_(INTEGER N,float[] A,INTEGER LDA,int[] IPIV,int[] JPIV,INTEGER INFO);
	public void sgetf2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void sgetrf_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void sgetrf2_(INTEGER M,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void sgetri_(INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgetrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void sggbak_(CHARACTER JOB,CHARACTER SIDE,INTEGER N,INTEGER ILO,INTEGER IHI,float[] LSCALE,float[] RSCALE,INTEGER M,float[] V,INTEGER LDV,INTEGER INFO);
	public void sggbal_(CHARACTER JOB,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER ILO,INTEGER IHI,float[] LSCALE,float[] RSCALE,float[] WORK,INTEGER INFO);
	public void sgges_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER SDIM,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] VSL,INTEGER LDVSL,float[] VSR,INTEGER LDVSR,float[] WORK,INTEGER LWORK,boolean BWORK,INTEGER INFO);
	public void sgges3_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER SDIM,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] VSL,INTEGER LDVSL,float[] VSR,INTEGER LDVSR,float[] WORK,INTEGER LWORK,boolean BWORK,INTEGER INFO);
	public void sggesx_(CHARACTER JOBVSL,CHARACTER JOBVSR,CHARACTER SORT,LOGICAL SELCTG,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER SDIM,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] VSL,INTEGER LDVSL,float[] VSR,INTEGER LDVSR,float[] RCONDE,float[] RCONDV,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,boolean BWORK,INTEGER INFO);
	public void sggev_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sggev3_(CHARACTER JOBVL,CHARACTER JOBVR,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sggevx_(CHARACTER BALANC,CHARACTER JOBVL,CHARACTER JOBVR,CHARACTER SENSE,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER ILO,INTEGER IHI,float[] LSCALE,float[] RSCALE,REAL ABNRM,REAL BBNRM,float[] RCONDE,float[] RCONDV,float[] WORK,INTEGER LWORK,int[] IWORK,boolean BWORK,INTEGER INFO);
	public void sggglm_(INTEGER N,INTEGER M,INTEGER P,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] D,float[] X,float[] Y,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgghd3_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgghrd_(CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,INTEGER INFO);
	public void sgglse_(INTEGER M,INTEGER N,INTEGER P,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] C,float[] D,float[] X,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sggqrf_(INTEGER N,INTEGER M,INTEGER P,float[] A,INTEGER LDA,float[] TAUA,float[] B,INTEGER LDB,float[] TAUB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sggrqf_(INTEGER M,INTEGER P,INTEGER N,float[] A,INTEGER LDA,float[] TAUA,float[] B,INTEGER LDB,float[] TAUB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sggsvd3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER N,INTEGER P,INTEGER K,INTEGER L,float[] A,INTEGER LDA,float[] B,INTEGER LDB,float[] ALPHA,float[] BETA,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] Q,INTEGER LDQ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void sggsvp3_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,float[] A,INTEGER LDA,float[] B,INTEGER LDB,REAL TOLA,REAL TOLB,INTEGER K,INTEGER L,float[] U,INTEGER LDU,float[] V,INTEGER LDV,float[] Q,INTEGER LDQ,int[] IWORK,float[] TAU,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgsvj0_(CHARACTER JOBV,INTEGER M,INTEGER N,float[] A,INTEGER LDA,float[] D,float[] SVA,INTEGER MV,float[] V,INTEGER LDV,REAL EPS,REAL SFMIN,REAL TOL,INTEGER NSWEEP,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgsvj1_(CHARACTER JOBV,INTEGER M,INTEGER N,INTEGER N1,float[] A,INTEGER LDA,float[] D,float[] SVA,INTEGER MV,float[] V,INTEGER LDV,REAL EPS,REAL SFMIN,REAL TOL,INTEGER NSWEEP,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sgtcon_(CHARACTER NORM,INTEGER N,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgtrfs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DLF,float[] DF,float[] DUF,float[] DU2,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgtsv_(INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] B,INTEGER LDB,INTEGER INFO);
	public void sgtsvx_(CHARACTER FACT,CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DLF,float[] DF,float[] DUF,float[] DU2,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sgttrf_(INTEGER N,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,INTEGER INFO);
	public void sgttrs_(CHARACTER TRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void sgtts2_(INTEGER ITRANS,INTEGER N,INTEGER NRHS,float[] DL,float[] D,float[] DU,float[] DU2,int[] IPIV,float[] B,INTEGER LDB);

}