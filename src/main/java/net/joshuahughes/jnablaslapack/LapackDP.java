package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDP extends Library
{

	public static LapackDP instance = (LapackDP) Native.loadLibrary("liblapack",LapackDP.class);

	public void dpbcon_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dpbequ_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void dpbrfs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dpbstf_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,INTEGER INFO);
	public void dpbsv_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,INTEGER INFO);
	public void dpbsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dpbtf2_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,INTEGER INFO);
	public void dpbtrf_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,INTEGER INFO);
	public void dpbtrs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,INTEGER INFO);
	public void dpftrf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,INTEGER INFO);
	public void dpftri_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,INTEGER INFO);
	public void dpftrs_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,double[] B,INTEGER LDB,INTEGER INFO);
	public void dpocon_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dpoequ_(INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void dpoequb_(INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void dporfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dporfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
	public void dposvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dposvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,int[] IWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE DLA_PORPVGRW);
	public void dpotf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dpotrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dpotrf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dpotri_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void dpotrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
	public void dppcon_(CHARACTER UPLO,INTEGER N,double[] AP,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dppequ_(CHARACTER UPLO,INTEGER N,double[] AP,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void dpprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dppsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,INTEGER INFO);
	public void dppsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dpptrf_(CHARACTER UPLO,INTEGER N,double[] AP,INTEGER INFO);
	public void dpptri_(CHARACTER UPLO,INTEGER N,double[] AP,INTEGER INFO);
	public void dpptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,INTEGER INFO);
	public void dpstf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] PIV,INTEGER RANK,DOUBLE TOL,double[] WORK,INTEGER INFO);
	public void dpstrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] PIV,INTEGER RANK,DOUBLE TOL,double[] WORK,INTEGER INFO);
	public void dptcon_(INTEGER N,double[] D,double[] E,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
	public void dpteqr_(CHARACTER COMPZ,INTEGER N,double[] D,double[] E,double[] Z,INTEGER LDZ,double[] WORK,INTEGER INFO);
	public void dptrfs_(INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] DF,double[] EF,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,INTEGER INFO);
	public void dptsv_(INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB,INTEGER INFO);
	public void dptsvx_(CHARACTER FACT,INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] DF,double[] EF,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,INTEGER INFO);
	public void dpttrf_(INTEGER N,double[] D,double[] E,INTEGER INFO);
	public void dpttrs_(INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB,INTEGER INFO);
	public void dptts2_(INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB);

}