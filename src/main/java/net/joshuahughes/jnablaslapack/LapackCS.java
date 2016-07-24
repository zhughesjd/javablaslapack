package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCS extends Library
{

	public static LapackCS instance = (LapackCS) Native.loadLibrary("liblapack",LapackCS.class);

	public void cspcon_(CHARACTER UPLO,INTEGER N,float[] AP,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,INTEGER INFO);
	public void cspmv_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] AP,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
	public void cspr_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] AP);
	public void csprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void cspsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void cspsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,float[] AFP,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void csptrf_(CHARACTER UPLO,INTEGER N,float[] AP,int[] IPIV,INTEGER INFO);
	public void csptri_(CHARACTER UPLO,INTEGER N,float[] AP,int[] IPIV,float[] WORK,INTEGER INFO);
	public void csptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] AP,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void csrscl_(INTEGER N,REAL SA,float[] SX,INTEGER INCX);
	public void cstedc_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void cstegr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,REAL ABSTOL,INTEGER M,float[] W,float[] Z,INTEGER LDZ,int[] ISUPPZ,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void cstein_(INTEGER N,float[] D,float[] E,INTEGER M,float[] W,int[] IBLOCK,int[] ISPLIT,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
	public void cstemr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,INTEGER M,float[] W,float[] Z,INTEGER LDZ,INTEGER NZC,int[] ISUPPZ,LOGICAL TRYRAC,float[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void csteqr_(CHARACTER COMPZ,INTEGER N,float[] D,float[] E,float[] Z,INTEGER LDZ,float[] WORK,INTEGER INFO);
	public void csycon_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,INTEGER INFO);
	public void csyconv_(CHARACTER UPLO,CHARACTER WAY,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] E,INTEGER INFO);
	public void csycon_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,REAL ANORM,REAL RCOND,float[] WORK,INTEGER INFO);
	public void csyequb_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,float[] S,REAL SCOND,REAL AMAX,float[] WORK,INTEGER INFO);
	public void csymv_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] A,INTEGER LDA,float[] X,INTEGER INCX,float[] BETA,float[] Y,INTEGER INCY);
	public void csyr_(CHARACTER UPLO,INTEGER N,float[] ALPHA,float[] X,INTEGER INCX,float[] A,INTEGER LDA);
	public void csyrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,float[] FERR,float[] BERR,float[] WORK,float[] RWORK,INTEGER INFO);
	public void csyrfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO);
	public void csysv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void csysvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,float[] FERR,float[] BERR,float[] WORK,INTEGER LWORK,float[] RWORK,INTEGER INFO);
	public void csysvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,float[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,float[] S,float[] B,INTEGER LDB,float[] X,INTEGER LDX,REAL RCOND,REAL RPVGRW,float[] BERR,INTEGER N_ERR_BNDS,float[] ERR_BNDS_NORM,float[] ERR_BNDS_COMP,INTEGER NPARAMS,float[] PARAMS,float[] WORK,float[] RWORK,INTEGER INFO,LOGICAL LSAME,REAL SLAMCH,REAL CLA_SYRPVGRW);
	public void csysv_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void csyswapr_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,INTEGER I1,INTEGER I2);
	public void csytf2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void csytf2_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void csytrf_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void csytrf_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void csytri_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER INFO);
	public void csytri2_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void csytri2x_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER NB,INTEGER INFO);
	public void csytri_rook_(CHARACTER UPLO,INTEGER N,float[] A,INTEGER LDA,int[] IPIV,float[] WORK,INTEGER INFO);
	public void csytrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);
	public void csytrs2_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,float[] WORK,INTEGER INFO);
	public void csytrs_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,float[] A,INTEGER LDA,int[] IPIV,float[] B,INTEGER LDB,INTEGER INFO);

}