package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZC extends Library
{

	public static LapackZC instance = (LapackZC) Native.loadLibrary("liblapack",LapackZC.class);

	public void zcgesv_(INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,int[] IPIV,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] WORK,float[] SWORK,double[] RWORK,INTEGER ITER,INTEGER INFO);
	public void zcposv_(CHARACTER UPLO,INTEGER N,INTEGER NRHS,double[] A,INTEGER LDA,double[] B,INTEGER LDB,double[] X,INTEGER LDX,double[] WORK,float[] SWORK,double[] RWORK,INTEGER ITER,INTEGER INFO);

}