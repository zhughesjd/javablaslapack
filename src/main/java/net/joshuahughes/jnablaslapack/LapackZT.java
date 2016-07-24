package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZT extends Library
{

	public static LapackZT instance = (LapackZT) Native.loadLibrary("liblapack",LapackZT.class);

	public void ztbcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,INTEGER KD,double[] AB,INTEGER LDAB,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztbrfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztbtrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER KD,INTEGER NRHS,double[] AB,INTEGER LDAB,double[] B,INTEGER LDB,INTEGER INFO);
	public void ztfsm_(CHARACTER TRANSR,CHARACTER SIDE,CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER M,INTEGER N,double[] ALPHA,double[] A,double[] B,INTEGER LDB);
	public void ztftri_(CHARACTER TRANSR,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER INFO);
	public void ztfttp_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] ARF,double[] AP,INTEGER INFO);
	public void ztfttr_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] ARF,double[] A,INTEGER LDA,INTEGER INFO);
	public void ztgevc_(CHARACTER SIDE,CHARACTER HOWMNY,boolean SELECT,INTEGER N,double[] S,INTEGER LDS,double[] P,INTEGER LDP,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztgex2_(LOGICAL WANTQ,LOGICAL WANTZ,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER J1,INTEGER INFO);
	public void ztgexc_(LOGICAL WANTQ,LOGICAL WANTZ,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER IFST,INTEGER ILST,INTEGER INFO);
	public void ztgsen_(INTEGER IJOB,LOGICAL WANTQ,LOGICAL WANTZ,boolean SELECT,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] ALPHA,double[] BETA,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,INTEGER M,DOUBLE PL,DOUBLE PR,double[] DIF,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER LIWORK,INTEGER INFO);
	public void ztgsja_(CHARACTER JOBU,CHARACTER JOBV,CHARACTER JOBQ,INTEGER M,INTEGER P,INTEGER N,INTEGER K,INTEGER L,double[] A,INTEGER LDA,double[] B,INTEGER LDB,DOUBLE TOLA,DOUBLE TOLB,double[] ALPHA,double[] BETA,double[] U,INTEGER LDU,double[] V,INTEGER LDV,double[] Q,INTEGER LDQ,double[] WORK,INTEGER NCYCLE,INTEGER INFO);
	public void ztgsna_(CHARACTER JOB,CHARACTER HOWMNY,boolean SELECT,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] S,double[] DIF,INTEGER MM,INTEGER M,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void ztgsy2_(CHARACTER TRANS,INTEGER IJOB,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,INTEGER LDC,double[] D,INTEGER LDD,double[] E,INTEGER LDE,double[] F,INTEGER LDF,DOUBLE SCALE,DOUBLE RDSUM,DOUBLE RDSCAL,INTEGER INFO);
	public void ztgsyl_(CHARACTER TRANS,INTEGER IJOB,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,INTEGER LDC,double[] D,INTEGER LDD,double[] E,INTEGER LDE,double[] F,INTEGER LDF,DOUBLE SCALE,DOUBLE DIF,double[] WORK,INTEGER LWORK,int[] IWORK,INTEGER INFO);
	public void ztpcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] AP,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztpmqrt_(CHARACTER SIDE,CHARACTER TRANS,INTEGER M,INTEGER N,INTEGER K,INTEGER L,INTEGER NB,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] WORK,INTEGER INFO);
	public void ztpqrt_(INTEGER M,INTEGER N,INTEGER L,INTEGER NB,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] T,INTEGER LDT,double[] WORK,INTEGER INFO);
	public void ztpqrt2_(INTEGER M,INTEGER N,INTEGER L,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] T,INTEGER LDT,INTEGER INFO);
	public void ztprfb_(CHARACTER SIDE,CHARACTER TRANS,CHARACTER DIRECT,CHARACTER STOREV,INTEGER M,INTEGER N,INTEGER K,INTEGER L,double[] V,INTEGER LDV,double[] T,INTEGER LDT,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] WORK,INTEGER LDWORK);
	public void ztprfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztptri_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] AP,INTEGER INFO);
	public void ztptrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] AP,double[] B,INTEGER LDB,INTEGER INFO);
	public void ztpttf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] AP,double[] ARF,INTEGER INFO);
	public void ztpttr_(CHARACTER UPLO,INTEGER N,double[] AP,double[] A,INTEGER LDA,INTEGER INFO);
	public void ztrcon_(CHARACTER NORM,CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,DOUBLE RCOND,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztrevc_(CHARACTER SIDE,CHARACTER HOWMNY,boolean SELECT,INTEGER N,double[] T,INTEGER LDT,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztrevc3_(CHARACTER SIDE,CHARACTER HOWMNY,boolean SELECT,INTEGER N,double[] T,INTEGER LDT,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,INTEGER LWORK,double[] RWORK,INTEGER LRWORK,INTEGER INFO);
	public void ztrexc_(CHARACTER COMPQ,INTEGER N,double[] T,INTEGER LDT,double[] Q,INTEGER LDQ,INTEGER IFST,INTEGER ILST,INTEGER INFO);
	public void ztrrfs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] FERR,double[] BERR,double[] WORK,double[] RWORK,INTEGER INFO);
	public void ztrsen_(CHARACTER JOB,CHARACTER COMPQ,boolean SELECT,INTEGER N,double[] T,INTEGER LDT,double[] Q,INTEGER LDQ,double[] W,INTEGER M,DOUBLE S,DOUBLE SEP,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void ztrsna_(CHARACTER JOB,CHARACTER HOWMNY,boolean SELECT,INTEGER N,double[] T,INTEGER LDT,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,double[] S,double[] SEP,INTEGER MM,INTEGER M,double[] WORK,INTEGER LDWORK,double[] RWORK,INTEGER INFO);
	public void ztrsyl_(CHARACTER TRANA,CHARACTER TRANB,INTEGER ISGN,INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] C,INTEGER LDC,DOUBLE SCALE,INTEGER INFO);
	public void ztrti2_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void ztrtri_(CHARACTER UPLO,CHARACTER DIAG,INTEGER N,double[] A,INTEGER LDA,INTEGER INFO);
	public void ztrtrs_(CHARACTER UPLO,CHARACTER TRANS,CHARACTER DIAG,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,INTEGER INFO);
	public void ztrttf_(CHARACTER TRANSR,CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] ARF,INTEGER INFO);
	public void ztrttp_(CHARACTER UPLO,INTEGER N,double[] A,INTEGER LDA,double[] AP,INTEGER INFO);
	public void ztzrzf_(INTEGER M,INTEGER N,double[] A,INTEGER LDA,double[] TAU,double[] WORK,INTEGER LWORK,INTEGER INFO);

}