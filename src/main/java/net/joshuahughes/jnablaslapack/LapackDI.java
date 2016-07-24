package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackDI extends Library
{

	public static LapackDI instance = (LapackDI) Native.loadLibrary("liblapack",LapackDI.class);

	public boolean disnan_(DOUBLE DIN);

}