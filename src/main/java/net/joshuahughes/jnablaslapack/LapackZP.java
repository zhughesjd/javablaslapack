package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZP extends Library
{

	public static LapackZP instance = (LapackZP) Native.loadLibrary("liblapack",LapackZP.class);

	public void zpbcon_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zpbequ_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void zpbrfs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zpbstf_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,INTEGER INFO);
	public void zpbsv_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,INTEGER INFO);
	public void zpbsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] AFB,INTEGER LDAFB,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zpbtf2_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,INTEGER INFO);
	public void zpbtrf_(CHARACTER UPLO,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,INTEGER INFO);
	public void zpbtrs_(CHARACTER UPLO,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,INTEGER INFO);
	public void zpftrf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,INTEGER INFO);
	public void zpftri_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,INTEGER INFO);
	public void zpftrs_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,double[] B,INTEGER LDB,INTEGER INFO);
	public void zpocon_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zpoequ_(INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void zpoequb_(INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void zporfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zporfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
	public void zposvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zposvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE ZLA_PORPVGRW);
	public void zpotf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void zpotrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void zpotrf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void zpotri_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void zpotrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
	public void zppcon_(CHARACTER UPLO,INTEGER N,double[] AP,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zppequ_(CHARACTER UPLO,INTEGER N,double[] AP,double[] S,DOUBLE SCOND,DOUBLE AMAX,INTEGER INFO);
	public void zpprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zppsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,INTEGER INFO);
	public void zppsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zpptrf_(CHARACTER UPLO,INTEGER N,double[] AP,INTEGER INFO);
	public void zpptri_(CHARACTER UPLO,INTEGER N,double[] AP,INTEGER INFO);
	public void zpptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,INTEGER INFO);
	public void zpstf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] PIV,INTEGER RANK,DOUBLE TOL,double[] WORK,INTEGER INFO);
	public void zpstrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] PIV,INTEGER RANK,DOUBLE TOL,double[] WORK,INTEGER INFO);
	public void zptcon_(INTEGER N,double[] D,double[] E,DOUBLE ANORM,DOUBLE RCOND,double[] RWORK,INTEGER INFO);
	public void zpteqr_(CHARACTER COMPZ,INTEGER N,double[] D,double[] E,double[] Z,INTEGER LDZ,double[] WORK,INTEGER INFO);
	public void zptrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] DF,double[] EF,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zptsv_(INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB,INTEGER INFO);
	public void zptsvx_(CHARACTER FACT,INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] DF,double[] EF,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zpttrf_(INTEGER N,double[] D,double[] E,INTEGER INFO);
	public void zpttrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB,INTEGER INFO);
	public void zptts2_(INTEGER IUPLO,INTEGER N,INTEGER NRHS,double[] D,double[] E,double[] B,INTEGER LDB);

}