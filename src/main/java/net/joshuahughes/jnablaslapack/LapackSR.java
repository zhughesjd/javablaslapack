package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSR extends Library
{

	public static LapackSR instance = (LapackSR) Native.loadLibrary("liblapack",LapackSR.class);

	public void srscl_(INTEGER N,REAL SA,float[] SX,INTEGER INCX);

}