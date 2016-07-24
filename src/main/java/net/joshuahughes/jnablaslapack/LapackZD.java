package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZD extends Library
{

	public static LapackZD instance = (LapackZD) Native.loadLibrary("liblapack",LapackZD.class);

	public void zdrscl_(INTEGER N,DOUBLE SA,double[] SX,INTEGER INCX);

}