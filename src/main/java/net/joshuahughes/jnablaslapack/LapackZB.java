package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZB extends Library
{

	public static LapackZB instance = (LapackZB) Native.loadLibrary("liblapack",LapackZB.class);

	public void zbbcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,INTEGER M,INTEGER P,INTEGER Q,double[] THETA,double[] PHI,double[] U1,INTEGER LDU1,double[] U2,INTEGER LDU2,double[] V1T,INTEGER LDV1T,double[] V2T,INTEGER LDV2T,double[] B11D,double[] B11E,double[] B12D,double[] B12E,double[] B21D,double[] B21E,double[] B22D,double[] B22E,double[] RWORK,INTEGER LRWORK,INTEGER INFO);
	public void zbdsqr_(CHARACTER UPLO,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,double[] D,double[] E,double[] VT,INTEGER LDVT,double[] U,INTEGER LDU,double[] C,INTEGER LDC,double[] RWORK,INTEGER INFO);

}