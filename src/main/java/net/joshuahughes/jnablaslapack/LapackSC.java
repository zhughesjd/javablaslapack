package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSC extends Library
{

	public static LapackSC instance = (LapackSC) Native.loadLibrary("liblapack",LapackSC.class);

	public float scsum1_(INTEGER N,float[] CX,INTEGER INCX);

}