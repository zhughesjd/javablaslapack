package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDR extends Library
{

	public static LapackDR instance = (LapackDR) Native.loadLibrary("liblapack",LapackDR.class);

	public void drscl_(INTEGER N,DOUBLE SA,double[] SX,INTEGER INCX);

}