package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCP extends Library
{

	public static LapackCP instance = (LapackCP) Native.loadLibrary("liblapack",LapackCP.class);

	public void cpbcon_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cpbequ_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void cpbrfs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cpbstf_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
	public void cpbsv_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
	public void cpbsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] AFB,INTEGER LDAFB,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cpbtf2_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
	public void cpbtrf_(CHARACTER UPLO,INTEGER N,INTEGER KD,float[] AB,INTEGER LDAB,INTEGER INFO);
	public void cpbtrs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,float[] AB,INTEGER LDAB,float[] B,INTEGER LDB,INTEGER INFO);
	public void cpftrf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER INFO);
	public void cpftri_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,float[] A,INTEGER INFO);
	public void cpftrs_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,float[] B,INTEGER LDB,INTEGER INFO);
	public void cpocon_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cpoequ_(INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void cpoequb_(INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void cporfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cporfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
	public void cposvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cposvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL CLA_PORPVGRW);
	public void cpotf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void cpotrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void cpotrf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void cpotri_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER INFO);
	public void cpotrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] B,INTEGER LDB,INTEGER INFO);
	public void cppcon_(CHARACTER UPLO,INTEGER N,float[] AP,REAL ANORM,REAL RCOND,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cppequ_(CHARACTER UPLO,INTEGER N,float[] AP,float[] S,REAL SCOND,REAL AMAX,INTEGER INFO);
	public void cpprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cppsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
	public void cppsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cpptrf_(CHARACTER UPLO,INTEGER N,float[] AP,INTEGER INFO);
	public void cpptri_(CHARACTER UPLO,INTEGER N,float[] AP,INTEGER INFO);
	public void cpptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] B,INTEGER LDB,INTEGER INFO);
	public void cpstf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] PIV,INTEGER RANK,REAL TOL,float[] WORK,INTEGER INFO);
	public void cpstrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] PIV,INTEGER RANK,REAL TOL,float[] WORK,INTEGER INFO);
	public void cptcon_(INTEGER N,float[] D,float[] E,REAL ANORM,REAL RCOND,float[] RWORK,INTEGER INFO);
	public void cpteqr_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
	public void cptrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] DF,float[] EF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cptsv_(INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,INTEGER INFO);
	public void cptsvx_(CHARACTER FACT,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] DF,float[] EF,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cpttrf_(INTEGER N,float[] D,float[] E,INTEGER INFO);
	public void cpttrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB,INTEGER INFO);
	public void cptts2_(INTEGER IUPLO,INTEGER N,INTEGER NRHS,float[] D,float[] E,float[] B,INTEGER LDB);

}