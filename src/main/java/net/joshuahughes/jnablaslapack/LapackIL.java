package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIL extends Library
{

	public static LapackIL instance = (LapackIL) Native.loadLibrary("liblapack",LapackIL.class);

	public int ilaclc_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
	public int ilaclr_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
	public int iladiag_(CHARACTER DIAG);
	public int iladlc_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);
	public int iladlr_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);
	public int ilaenv_(INTEGER ISPEC,CHARACTER NAME,CHARACTER OPTS,INTEGER N1,INTEGER N2,INTEGER N3,INTEGER N4);
	public int ilaprec_(CHARACTER PREC);
	public int ilaslc_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
	public int ilaslr_(INTEGER M,INTEGER N,float[] A,INTEGER LDA);
	public int ilatrans_(CHARACTER TRANS);
	public int ilauplo_(CHARACTER UPLO);
	public void ilaver_(INTEGER VERS_MAJOR,INTEGER VERS_MINOR,INTEGER VERS_PATCH);
	public int ilazlc_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);
	public int ilazlr_(INTEGER M,INTEGER N,double[] A,INTEGER LDA);

}