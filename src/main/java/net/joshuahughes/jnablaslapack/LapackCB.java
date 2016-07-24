package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCB extends Library
{

	public static LapackCB instance = (LapackCB) Native.loadLibrary("liblapack",LapackCB.class);

	public void cbbcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,INTEGER M,INTEGER P,INTEGER Q,float[] THETA,float[] PHI,float[] U1,INTEGER LDU1,float[] U2,INTEGER LDU2,float[] V1T,INTEGER LDV1T,float[] V2T,INTEGER LDV2T,float[] B11D,float[] B11E,float[] B12D,float[] B12E,float[] B21D,float[] B21E,float[] B22D,float[] B22E,float[] RWORK,INTEGER LRWORK,INTEGER INFO);
	public void cbdsqr_(CHARACTER UPLO,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,float[] D,float[] E,float[] VT,INTEGER LDVT,float[] U,INTEGER LDU,float[] C,INTEGER LDC,float[] RWORK,INTEGER INFO);

}