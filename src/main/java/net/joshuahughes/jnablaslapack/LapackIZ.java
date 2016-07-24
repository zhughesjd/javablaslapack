package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIZ extends Library
{

	public static LapackIZ instance = (LapackIZ) Native.loadLibrary("liblapack",LapackIZ.class);

	public int izmax1_(INTEGER N,double[] ZX,INTEGER INCX);

}