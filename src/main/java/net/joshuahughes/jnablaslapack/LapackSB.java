package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSB extends Library
{

	public static LapackSB instance = (LapackSB) Native.loadLibrary("liblapack",LapackSB.class);

	public void sbbcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,INTEGER M,INTEGER P,INTEGER Q,float[] THETA,float[] PHI,float[] U1,INTEGER LDU1,float[] U2,INTEGER LDU2,float[] V1T,INTEGER LDV1T,float[] V2T,INTEGER LDV2T,float[] B11D,float[] B11E,float[] B12D,float[] B12E,float[] B21D,float[] B21E,float[] B22D,float[] B22E,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void sbdsdc_(CHARACTER UPLO,CHARACTER COMPQ,INTEGER N,float[] D,float[] E,float[] U,INTEGER LDU,float[] VT,INTEGER LDVT,float[] Q,int[] IQ,float[] WORK,int[] IWORK,INTEGER INFO);
	public void sbdsqr_(CHARACTER UPLO,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,float[] D,float[] E,float[] VT,INTEGER LDVT,float[] U,INTEGER LDU,float[] C,INTEGER LDC,float[] WORK,INTEGER INFO);
	public void sbdsvdx_(CHARACTER UPLO,CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,float[] D,float[] E,REAL VL,REAL VU,INTEGER IL,INTEGER IU,INTEGER NS,float[] S,float[] Z,INTEGER LDZ,float[] WORK,int[] IWORK,INTEGER INFO);

}