package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackZR extends Library
{

	public static LapackZR instance = (LapackZR) Native.loadLibrary("liblapack",LapackZR.class);

	public void zrot_(INTEGER N,double[] CX,INTEGER INCX,double[] CY,INTEGER INCY,DOUBLE C,double[] S);

}