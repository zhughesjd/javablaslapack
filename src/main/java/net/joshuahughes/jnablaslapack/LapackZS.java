package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZS extends Library
{

	public static LapackZS instance = (LapackZS) Native.loadLibrary("liblapack",LapackZS.class);

	public void zspcon_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
	public void zspmv_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] AP,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
	public void zspr_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] AP);
	public void zsprfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zspsv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void zspsvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,double[] AFP,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zsptrf_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,INTEGER INFO);
	public void zsptri_(CHARACTER UPLO,INTEGER N,double[] AP,int[] IPIV,double[] WORK,INTEGER INFO);
	public void zsptrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] AP,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void zstedc_(CHARACTER COMPZ,INTEGER N,double[] D,double[] E,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void zstegr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,double[] D,double[] E,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,DOUBLE ABSTOL,INTEGER M,double[] W,double[] Z,INTEGER LDZ,int[] ISUPPZ,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void zstein_(INTEGER N,double[] D,double[] E,INTEGER M,double[] W,int[] IBLOCK,int[] ISPLIT,double[] Z,INTEGER LDZ,double[] WORK,int[] IWORK,int[] IFAIL,INTEGER INFO);
	public void zstemr_(CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,double[] D,double[] E,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,INTEGER M,double[] W,double[] Z,INTEGER LDZ,INTEGER NZC,int[] ISUPPZ,LOGICAL TRYRAC,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void zsteqr_(CHARACTER COMPZ,INTEGER N,double[] D,double[] E,double[] Z,INTEGER LDZ,double[] WORK,INTEGER INFO);
	public void zsycon_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
	public void zsyconv_(CHARACTER UPLO,CHARACTER WAY,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] E,INTEGER INFO);
	public void zsycon_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,DOUBLE ANORM,DOUBLE RCOND,double[] WORK,INTEGER INFO);
	public void zsyequb_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] S,DOUBLE SCOND,DOUBLE AMAX,double[] WORK,INTEGER INFO);
	public void zsymv_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] A,INTEGER LDA,double[] X,INTEGER INCX,double[] BETA,double[] Y,INTEGER INCY);
	public void zsyr_(CHARACTER UPLO,INTEGER N,double[] ALPHA,double[] X,INTEGER INCX,double[] A,INTEGER LDA);
	public void zsyrfs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zsyrfsx_(CHARACTER UPLO,CHARACTER EQUED,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO);
	public void zsysv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void zsysvx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,double[] FERR,double[] BERR,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER INFO);
	public void zsysvxx_(CHARACTER FACT,CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] AF,INTEGER LDAF,int[] IPIV,CHARACTER EQUED,double[] S,double[] B,INTEGER LDB,double[] X,INTEGER LDX,DOUBLE RCOND,DOUBLE RPVGRW,double[] BERR,INTEGER N_ERR_BNDS,double[] ERR_BNDS_NORM,double[] ERR_BNDS_COMP,INTEGER NPARAMS,double[] PARAMS,double[] WORK,double[] RWORK,INTEGER INFO,LOGICAL LSAME,DOUBLE DLAMCH,DOUBLE ZLA_SYRPVGRW);
	public void zsysv_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void zsyswapr_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,INTEGER I1,INTEGER I2);
	public void zsytf2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void zsytf2_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,INTEGER INFO);
	public void zsytrf_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void zsytrf_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void zsytri_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER INFO);
	public void zsytri2_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void zsytri2x_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER NB,INTEGER INFO);
	public void zsytri_rook_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,int[] IPIV,double[] WORK,INTEGER INFO);
	public void zsytrs_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);
	public void zsytrs2_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] WORK,INTEGER INFO);
	public void zsytrs_rook_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,INTEGER INFO);

}