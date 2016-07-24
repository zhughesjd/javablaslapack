package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackIE extends Library
{

	public static LapackIE instance = (LapackIE) Native.loadLibrary("liblapack",LapackIE.class);

	public int ieeeck_(INTEGER ISPEC,REAL ZERO,REAL ONE);

}