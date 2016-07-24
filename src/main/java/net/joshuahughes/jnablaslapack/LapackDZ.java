package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDZ extends Library
{

	public static LapackDZ instance = (LapackDZ) Native.loadLibrary("liblapack",LapackDZ.class);

	public double dzsum1_(INTEGER N,double[] CX,INTEGER INCX);

}