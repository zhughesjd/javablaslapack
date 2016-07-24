package net.joshuahughes.jnablaslapack;

import com.sun.jna.Library;
import com.sun.jna.Native;
import net.joshuahughes.jnablaslapack.pointer.*;

public interface LapackLS extends Library
{

	public static LapackLS instance = (LapackLS) Native.loadLibrary("liblapack",LapackLS.class);

	public boolean lsamen_(INTEGER N,CHARACTER CA,CHARACTER CB);

}