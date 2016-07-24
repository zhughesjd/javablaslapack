package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackSI extends Library
{

	public static LapackSI instance = (LapackSI) Native.loadLibrary("liblapack",LapackSI.class);

	public boolean sisnan_(REAL SIN);

}