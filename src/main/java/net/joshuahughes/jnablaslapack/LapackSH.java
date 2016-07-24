package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSH extends Library
{

	public static LapackSH instance = (LapackSH) Native.loadLibrary("liblapack",LapackSH.class);

	public void shgeqz_(CHARACTER JOB,CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] T,INTEGER LDT,float[] ALPHAR,float[] ALPHAI,float[] BETA,float[] Q,INTEGER LDQ,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);
	public void shsein_(CHARACTER SIDE,CHARACTER EIGSRC,CHARACTER INITV,boolean SELECT,INTEGER N,float[] H,INTEGER LDH,float[] WR,float[] WI,float[] VL,INTEGER LDVL,float[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,float[] WORK,int[] IFAILL,int[] IFAILR,INTEGER INFO);
	public void shseqr_(CHARACTER JOB,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,float[] H,INTEGER LDH,float[] WR,float[] WI,float[] Z,INTEGER LDZ,float[] WORK,INTEGER LWORK,INTEGER INFO);

}