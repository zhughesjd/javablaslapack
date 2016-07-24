package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSP extends Library
{

	public static LapackSP instance = (LapackSP) Native.loadLibrary("liblapack",LapackSP.class);

	public void spbcon_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
	public void spbequ_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void spbrfs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void spbstf_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
	public void spbsv_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
	public void spbsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void spbtf2_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
	public void spbtrf_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
	public void spbtrs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
	public void spftrf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER INFO);
	public void spftri_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER INFO);
	public void spftrs_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,float[] B,INTEGER LDB,INTEGER INFO);
	public void spocon_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
	public void spoequ_(INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void spoequb_(INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void sporfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sporfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
	public void sposvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sposvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL SLA_PORPVGRW);
	public void spotf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void spotrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void spotrf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void spotri_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void spotrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
	public void sppcon_(CHARACTER UPLO,INTEGER N,float[] AP,REAL ANORM,REAL RCOND,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sppequ_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void spprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sppsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
	public void sppsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,int[] IWORK,INTEGER INFO);
	public void spptrf_(CHARACTER UPLO,INTEGER N,float[] AP,INTEGER INFO);
	public void spptri_(CHARACTER UPLO,INTEGER N,float[] AP,INTEGER INFO);
	public void spptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
	public void spstf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] PIV,INTEGER RANK,REAL TOL,float[] WORK,INTEGER INFO);
	public void spstrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] PIV,INTEGER RANK,REAL TOL,float[] WORK,INTEGER INFO);
	public void sptcon_(INTEGER N,float[] D,float[] E,REAL ANORM,REAL RCOND,float[] WORK,INTEGER INFO);
	public void spteqr_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
	public void sptrfs_(INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] DF,float[] EF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,INTEGER INFO);
	public void sptsv_(INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,INTEGER INFO);
	public void sptsvx_(CHARACTER FACT,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] DF,float[] EF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,INTEGER INFO);
	public void spttrf_(INTEGER N,float[] D,float[] E,INTEGER INFO);
	public void spttrs_(INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,INTEGER INFO);
	public void sptts2_(INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB);

}