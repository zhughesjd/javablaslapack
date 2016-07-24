package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDB extends Library
{

	public static LapackDB instance = (LapackDB) Native.loadLibrary("liblapack",LapackDB.class);

	public void dbbcsd_(CHARACTER JOBU1,CHARACTER JOBU2,CHARACTER JOBV1T,CHARACTER JOBV2T,CHARACTER TRANS,INTEGER M,INTEGER P,INTEGER Q,double[] THETA,double[] PHI,double[] U1,INTEGER LDU1,double[] U2,INTEGER LDU2,double[] V1T,INTEGER LDV1T,double[] V2T,INTEGER LDV2T,double[] B11D,double[] B11E,double[] B12D,double[] B12E,double[] B21D,double[] B21E,double[] B22D,double[] B22E,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dbdsdc_(CHARACTER UPLO,CHARACTER COMPQ,INTEGER N,double[] D,double[] E,double[] U,INTEGER LDU,double[] VT,INTEGER LDVT,double[] Q,int[] IQ,double[] WORK,int[] IWORK,INTEGER INFO);
	public void dbdsqr_(CHARACTER UPLO,INTEGER N,INTEGER NCVT,INTEGER NRU,INTEGER NCC,double[] D,double[] E,double[] VT,INTEGER LDVT,double[] U,INTEGER LDU,double[] C,INTEGER LDC,double[] WORK,INTEGER INFO);
	public void dbdsvdx_(CHARACTER UPLO,CHARACTER JOBZ,CHARACTER RANGE,INTEGER N,double[] D,double[] E,DOUBLE VL,DOUBLE VU,INTEGER IL,INTEGER IU,INTEGER NS,double[] S,double[] Z,INTEGER LDZ,double[] WORK,int[] IWORK,INTEGER INFO);

}