package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIC extends Library
{

	public static LapackIC instance = (LapackIC) Native.loadLibrary("liblapack",LapackIC.class);

	public int icmax1_(INTEGER N,float[] CX,INTEGER INCX);

}