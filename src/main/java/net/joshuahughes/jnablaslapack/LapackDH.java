package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDH extends Library
{

	public static LapackDH instance = (LapackDH) Native.loadLibrary("liblapack",LapackDH.class);

	public void dhgeqz_(CHARACTER JOB,CHARACTER COMPQ,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] T,INTEGER LDT,double[] ALPHAR,double[] ALPHAI,double[] BETA,double[] Q,INTEGER LDQ,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);
	public void dhsein_(CHARACTER SIDE,CHARACTER EIGSRC,CHARACTER INITV,boolean SELECT,INTEGER N,double[] H,INTEGER LDH,double[] WR,double[] WI,double[] VL,INTEGER LDVL,double[] VR,INTEGER LDVR,INTEGER MM,INTEGER M,double[] WORK,int[] IFAILL,int[] IFAILR,INTEGER INFO);
	public void dhseqr_(CHARACTER JOB,CHARACTER COMPZ,INTEGER N,INTEGER ILO,INTEGER IHI,double[] H,INTEGER LDH,double[] WR,double[] WI,double[] Z,INTEGER LDZ,double[] WORK,INTEGER LWORK,INTEGER INFO);

}