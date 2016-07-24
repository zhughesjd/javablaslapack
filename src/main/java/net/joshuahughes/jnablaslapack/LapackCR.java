package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackCR extends Library
{

	public static LapackCR instance = (LapackCR) Native.loadLibrary("liblapack",LapackCR.class);

	public void crot_(INTEGER N,float[] CX,INTEGER INCX,float[] CY,INTEGER INCY,REAL C,float[] S);

}